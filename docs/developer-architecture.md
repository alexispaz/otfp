# Developer Architecture

This page documents the current architecture before refactoring. The refactor
branch does not need to preserve backward compatibility, but it should preserve
scientific behavior unless a feature is explicitly removed.

## Current Layers

| Layer | Files | Responsibility |
| --- | --- | --- |
| NAMD callback | `calcforces.tcl` | Register NAMD groups, load coordinates, call C/Tcl update flow, return forces and energy. |
| Tcl setup | `cfacv.tcl`, `cvs.tcl`, `restraints.tcl`, `chapeau.tcl`, `data.tcl` | Parse user variables, allocate C objects through SWIG, and prepare runtime state. |
| SWIG API | `cfacv.i`, `cvs.i`, public headers | Expose C structs/functions to Tcl. |
| Runtime public API | `cfacv.h` | Expose opaque runtime, restraint, and optimizer types plus user/SWIG-facing setup and update functions. |
| Runtime internals | `cfacv_internal.h`, `cfacv.c` | Define `DataSpace`, `restraint`, optimizer structs, and orchestrate CVs, forces, targets, and FES updates. |
| DataSpace ownership | `cfacv_dataspace.c` | Allocate/free `DataSpace`, register atom centers, CVs, restraints, and chapeau objects. |
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
| Tools | `catbinsp.c`, `simetrize.c`, `convchapeau.c`, `addchapeaus.c`, `cropchapeaus.c` | Offline manipulation of chapeau outputs/states. |

## Main Runtime Flow

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
- `cvs.h` declares `typedef struct cv cv` plus public CV helpers.
- `cvs_internal.h` defines `struct cv` and the internal calculation/boundary
  functions used by the C runtime.
- `cfacv.h` declares opaque `DataSpace`, `restraint`, and optimizer types plus
  the Tcl-facing runtime setup/update API.
- `cfacv_internal.h` defines the runtime structs and internal force/evolution
  helpers used by C modules.

The remaining public headers are intentionally narrower than the old `_obj`
headers, but the Tcl/SWIG surface still includes more functions than a final
user API should expose:

- `cfacv.i` includes `cfacv.h`.
- `cvs.i` includes the public `chapeau.h` and public `cvs.h`.

During refactor, split headers into:

- Public Tcl/SWIG API headers.
- Internal headers for C modules.
- Tool-only headers if needed.

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
| `chapeau_internal.h` | Internal chapeau state is now isolated but still shared broadly among C modules and tools. | Replace direct field access with focused helper APIs where it reduces coupling. |
| `chapeau_sparse.c` | Duplicate implementation marked under construction. | Remove from default architecture or replace with a backend interface. |
| SWIG headers | Internal helpers become Tcl-visible API. | Expose a small deliberate API. |
| Error handling | Most failures call `exit`, which is hard to test and hard to recover from Tcl. | Return status codes in C and convert to Tcl errors at the boundary. |
| Allocation | Constructors do not consistently initialize optional pointers; most types lack frees. | Introduce zero-initializing constructors and destructors. |

## Refactor Sequence

1. Narrow the Tcl/SWIG API to deliberate user-facing functions.
2. Replace broad direct struct access with focused helper APIs where it reduces coupling.
3. Decide explicitly which experimental/legacy features remain.

## Behavior Guardrails

Fast tests should pass after each mechanical refactor:

```sh
meson compile -C build
meson install -C build
meson test -C build examples
```

Documentation should be updated in the same commit when a user-visible feature,
configuration key, or developer boundary changes.
