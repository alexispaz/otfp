# Developer Architecture

This page documents the current architecture during refactoring. The refactor
branch does not need to preserve backward compatibility, but it should preserve
scientific behavior unless a feature is explicitly removed.

## Current Layers

| Layer | Files | Responsibility |
| --- | --- | --- |
| NAMD callback | `calcforces.tcl` | Register NAMD groups, load coordinates, call C/Tcl update flow, return forces and energy. |
| Tcl setup | `cfacv.tcl`, `cvs.tcl`, `restraints.tcl`, `chapeau.tcl`, `data.tcl` | Parse user variables, allocate C objects through SWIG, and prepare runtime state. |
| SWIG API | `cfacv.i`, `cvs.i`, `cfacv_tcl.h`, `cvs_tcl.h` | Expose a deliberate Tcl-facing C API. |
| Runtime public API | `cfacv.h` | Expose opaque runtime, restraint, and optimizer types plus user/SWIG-facing setup and update functions. |
| Runtime internals | `cfacv_internal.h`, `cfacv.c` | Define `DataSpace`, `restraint`, optimizer structs, and expose runtime helpers for CV computation, energy, and reporting. `DataSpace` fields use explicit total/count/ownership names for core readability. |
| Runtime force loop | `cfacv_force.c` | Orchestrate per-step force clearing, restraint values, restraint forces, CV boundary forces, FES/chapeau updates, and target evolution. |
| Driver embedding facade | `otfp_engine.h`, `otfp_engine.c` | Provide a small C API for non-Tcl drivers such as GEMS. |
| DataSpace ownership | `cfacv_dataspace.c` | Allocate/free `DataSpace`, register atom centers, CVs, restraints, and FES/chapeau objects. |
| Restraints | `cfacv_restraint.c` | Restraint construction, restart state, energy functions, boundary functions, and TAMD/SMD/ADAM target evolution. |
| CV public API | `cvs.h` | Expose the opaque `cv` type and user/SWIG-facing CV lifecycle, type, reference, and setup helpers. |
| CV internals | `cvs_internal.h`, `cvs.c`, `measurements.c` | Define the concrete `cv` layout, CV descriptors, boundary functions, and value/gradient calculations. |
| FES/chapeau public API | `chapeau.h` | Expose the opaque `chapeau` type and user/SWIG-facing lifecycle, output, and restart functions. |
| FES/chapeau internals | `chapeau_internal.h`, `chapeau.c` | Define the concrete `chapeau` layout and internal module contracts for lifecycle, grid metadata, and size/grid comparisons. |
| FES/chapeau grid utilities | `chapeau_grid.c` | Crop, sum, and symmetrize chapeau grid statistics. |
| FES/chapeau accumulation | `chapeau_accumulate.c` | Accumulate 1D/2D force statistics into chapeau matrices and hit masks. |
| FES/chapeau IO | `chapeau_io.c` | Manage binary output, restart state, lambda loading, and replica serialization. |
| FES/chapeau solver | `chapeau_solver.c` | Build reduced band matrices and solve chapeau coefficients with LAPACK. |
| FES/chapeau evaluation | `chapeau_eval.c` | Evaluate 1D/2D simplex interpolation, gradients, and peak-file loading. |
| Tools | `catbinsp.c`, `simetrize.c`, `addchapeaus.c`, `cropchapeaus.c` | Offline manipulation of chapeau outputs/states. |

## Driver Boundaries

NAMD is currently the only complete production driver in this repository, but
it is not the engine model. The desired architecture is:

| Boundary | Responsibility |
| --- | --- |
| Core engine | Own `DataSpace`, CVs, restraints, FES/chapeau state, per-step force/energy calculation, and restart/output scheduling. |
| NAMD/Tcl adapter | Translate NAMD TclForces callbacks and user Tcl variables into core engine calls. |
| Mock adapter | Emulate enough of NAMD TclForces for fast regression tests. |
| Future GEMS adapter | Bind the `otfp_engine` C facade from Fortran with `iso_c_binding`, without going through Tcl or SWIG. |

The GEMS integration plan is documented in
[GEMS integration](integration-gems.md). Refactors should keep moving NAMD/Tcl
assumptions out of the core engine boundary.

## Current NAMD Runtime Flow

```text
NAMD calcforces
  -> loadcoords p
  -> Tcl_UpdateDataSpace
     -> DataSpace_centerPos copies positions into DataSpace.R
     -> DataSpace_ComputeCVs
        -> cv.calc for every CV
     -> DataSpace_RestrainingForces
        -> compute restraint values from CV linear combinations
        -> compute restraint energy/force
        -> accumulate center forces into DataSpace.R
        -> apply CV boundary forces
        -> accumulate and solve chapeau state when FES is enabled
        -> evolve TAMD/SMD/ADAM restraint targets
        -> save restart state
     -> addforce for each NAMD group
     -> addenergy from DataSpace_RestraintEnergy
```

## Public API Boundaries

The chapeau API has been split so SWIG sees only the opaque public interface:

- `chapeau.h` declares `typedef struct CHAPEAU chapeau` plus public lifecycle,
  output, and restart functions.
- `chapeau_internal.h` defines `struct CHAPEAU` and the internal module
  functions used by C modules and command-line tools.
- Runtime orchestration uses chapeau helper functions for sampling,
  accumulation, solve scheduling, output scheduling, and default restart saves
  instead of writing chapeau fields directly.
- `cvs.h` declares `typedef struct cv cv` plus public CV helpers.
- `cvs_internal.h` defines `struct cv`, type lookup helpers, and the internal
  calculation/boundary functions used by the C runtime.
- Runtime force orchestration uses CV helper functions for compute, value,
  gradient, boundary force, and output scheduling; direct CV field access is
  limited to CV implementation code.
- `cfacv.h` declares opaque `DataSpace`, `restraint`, and optimizer types plus
  the supported C runtime setup/update API.
- `cfacv_internal.h` defines the runtime structs and internal force/evolution
  helpers used by C modules. It also owns low-level constructors, type lookup
  helpers, debug dump helpers, and file wrappers that are not part of the public
  C or Tcl API.
- Runtime force orchestration uses restraint helper functions for value
  accumulation, energy/force calculation, output scheduling, chapeau sampling,
  target evolution, and energy reads; direct restraint field access is limited
  to restraint implementation and restart serialization.

The SWIG surface no longer includes the broad C headers directly:

- `cfacv.i` includes `cfacv_tcl.h`.
- `cvs.i` includes `cvs_tcl.h`.

These Tcl headers declare only the functions used by the main Tcl setup/runtime
flow: `DataSpace` lifecycle and setup, CV setup helpers, restraint setup and
target updates, per-step compute/force/energy calls, and chapeau restart calls.
They intentionally do not export lower-level constructors, type lookup helpers,
debug dump helpers, or internal chapeau evaluation/serialization routines.

The remaining broad C headers are still used for C compilation and tests.
During the next refactor passes, continue splitting them into:

- Internal headers for C modules.
- Tool-only headers if needed.
- Narrow C headers for tests and external embedding, if that remains a supported
  use case.

This is the first step before deleting or renaming functions.

## Data Ownership

Current allocated objects:

- `DataSpace`
- `atomCenterStruct`
- `cv`
- `restraint`
- `tamdOptStruct`
- `smdOptStruct`
- `adamOptStruct`
- `chapeau`

The refactor branch now has destructors for the core runtime objects:

- `atomCenterStruct_free`
- `centerStruct_free`
- `cv_free`
- `restraint_free`
- `DataSpace_free`
- `chapeau_free`

The `lifecycle` Meson test exercises the basic constructor/destructor path for
`DataSpace`, atom centers, one CV, and one restraint. Broader ownership tests
should be added as modules are split.

## High-Value Refactor Targets

| Target | Problem | Suggested direction |
| --- | --- | --- |
| `cfacv.c` | Mixes construction, force loop, output, restart, target evolution, and OTFP solve orchestration. | Split into `dataspace`, `restraint`, `evolution`, and `force_loop` modules. |
| `cvs.c` | Switch-heavy factory, global setup state for ZSD CVs, and all CV implementations in one file. | Use CV descriptor tables and per-CV parameter structs. |
| Internal state headers | Chapeau, CV, and restraint state are isolated, and runtime orchestration now uses focused helper APIs for their hot paths and diagnostics. Tools, serialization, and implementation modules still access fields directly. | Keep narrowing direct field access where it reduces coupling without obscuring numerical code. |
| SWIG headers | Internal helpers become Tcl-visible API. | Expose a small deliberate API. |
| Error handling | Most failures call `exit`, which is hard to test and hard to recover from Tcl. | Return status codes in C and convert to Tcl errors at the boundary. |
| Allocation | Constructors do not consistently initialize optional pointers; most types lack frees. | Introduce zero-initializing constructors and destructors. |

## Refactor Sequence

1. Narrow the Tcl/SWIG API to deliberate user-facing functions.
2. Replace broad direct struct access with focused helper APIs where it reduces coupling.
3. Expand the driver-oriented C embedding facade for GEMS and non-NAMD tests.
4. Decide explicitly which experimental/legacy features remain.

## Behavior Guardrails

Fast tests should pass after each mechanical refactor:

```sh
meson compile -C build
meson install -C build
meson test -C build examples
```

Documentation should be updated in the same commit when a user-visible feature,
configuration key, or developer boundary changes.
