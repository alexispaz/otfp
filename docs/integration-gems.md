# GEMS Integration Plan

This page describes how to prepare OTFP for integration into
[GEMS](https://github.com/alexispaz/GEMS) as an in-tree free-energy engine.

The goal is not to preserve OTFP as a standalone community dependency. OTFP is
stable laboratory code, while GEMS is the active experimentation environment for
free-energy methods and non-NAMD potentials such as DFTB and tight-binding
models. The practical integration path is to vendor the OTFP core into GEMS and
drive it from GEMS, while keeping the NAMD/Tcl path as an adapter.

## Current Shape

OTFP currently has four conceptual parts:

| Part | Current files | Role |
| --- | --- | --- |
| Core engine | `cfacv*.c`, `cvs*.c`, `chapeau*.c`, public/internal headers | Owns CVs, restraints, chapeau accumulation, solve scheduling, restart state, and force/energy calculation. |
| NAMD adapter | `calcforces.tcl`, `cfacv.tcl`, `cvs.tcl`, `restraints.tcl`, `chapeau.tcl`, `data.tcl` | Converts NAMD TclForces data and user Tcl variables into core engine calls. |
| Tcl/SWIG boundary | `cfacv.i`, `cvs.i`, `cfacv_tcl.h`, `cvs_tcl.h` | Exposes the selected C API to Tcl. This is an adapter boundary, not the engine model. |
| Mock adapter | `mock.tcl`, `examples/test.sh`, selected examples | Emulates enough NAMD TclForces behavior to run fast regression tests without a full NAMD calculation. |

For GEMS, the core engine is the useful part. The Tcl/SWIG and NAMD layers are
reference adapters and regression tools.

## Engine Model

A molecular driver does not need to know about Tcl or NAMD. It needs a
step-level OTFP engine with this lifecycle:

1. Allocate one `DataSpace`.
2. Register atom centers with masses and atom indices.
3. Register CVs.
4. Register restraints.
5. Optionally register chapeau/FES state.
6. On each step, copy current center positions into `DataSpace`.
7. Compute CV values.
8. Compute restraint forces, chapeau accumulation, target evolution, and
   restart/output scheduling.
9. Read forces from `DataSpace` center buffers and add them to the driver.
10. Read restraint energy and add it to the driver energy accounting.
11. Free the `DataSpace` and owned objects at shutdown.

The current public C calls that match this model are:

| Stage | Current C API |
| --- | --- |
| Allocate/free | `NewDataSpace`, `DataSpace_free` |
| Register centers | `DataSpace_AddAtomCenter`, `DataSpace_centerPos` |
| Register CVs | `DataSpace_add_cv`, `set_zsd_circle`, `set_zsd_ring`, `set_d2chap`, `set_line`, `cv_access_ref`, `cv_access_ref2` |
| Register restraints | `DataSpace_AddRestr`, `restr_AddTamdOpt`, `restr_AddSmdOpt`, `restr_AddAdamOpt` |
| Register FES | `DataSpace_SetupChapeau`, `DataSpace_get_chapeauadress` |
| Step update | `DataSpace_ComputeCVs`, `DataSpace_RestrainingForces` |
| Read results | `DataSpace_centerPos`, `DataSpace_RestraintEnergy`, `restr_getz`, `restr_getu` |
| Restart hooks | `chapeau_savestate`, `chapeau_loadstate`, `chapeau_loadlambda` |

This API is usable, but it still reflects the historical Tcl/NAMD workflow.
Before GEMS integration, prefer adding a small embedding facade instead of
teaching GEMS every historical function.

## Recommended GEMS Strategy

Use OTFP as vendored source inside GEMS, not as an external installed
dependency:

- GEMS is the active laboratory for the methods.
- OTFP has limited standalone maintenance pressure.
- Keeping the code in-tree makes it easier to couple OTFP with GEMS potentials,
  integrators, restart formats, and future free-energy experiments.
- The NAMD adapter can remain as an optional compatibility path or regression
  reference.

Do not translate the whole OTFP codebase to Fortran before integration. The
lower-risk path is:

1. Keep the numerical core in C.
2. Add a small C embedding facade with stable names and ownership rules.
3. Bind that facade from Fortran using `iso_c_binding`.
4. Translate or rewrite pieces later only when the integration shows a clear
   benefit.

This avoids mixing language conversion, build-system migration, memory
ownership changes, and scientific validation in one step.

## Proposed Embedding Facade

Add a driver-oriented API after the current headers are stable enough. A future
facade could live in files such as `otfp_engine.h` and `otfp_engine.c`.

The facade should hide Tcl/NAMD naming and expose a small workflow:

```c
typedef struct otfp_engine otfp_engine;

otfp_engine *otfp_engine_create(int ncenters, int ncvs, int nrestr, long seed);
void otfp_engine_free(otfp_engine *engine);

int otfp_engine_add_center(otfp_engine *engine, int natoms, const int *indices,
                           const double *masses);
int otfp_engine_set_center_position(otfp_engine *engine, int center,
                                    const double xyz[3]);
const double *otfp_engine_center_force(otfp_engine *engine, int center);

int otfp_engine_compute(otfp_engine *engine, int first_step, int step,
                        double external_bias);
double otfp_engine_bias_energy(otfp_engine *engine);
```

That shape is easier to bind from Fortran than the historical Tcl-facing API.
It also gives GEMS one stable dependency point while the internals continue to
be refactored.

## Fortran Binding Notes

The first GEMS integration should use `iso_c_binding` against the embedding
facade:

- Store the opaque C handle as `type(c_ptr)`.
- Pass scalar integers as `integer(c_int)` or the exact width chosen by the C
  facade.
- Pass coordinates and forces as contiguous `real(c_double)` arrays.
- Keep allocation and destruction on the C side.
- Convert Fortran strings to C strings only at setup time; avoid string-heavy
  calls in the per-step path.
- Document index base rules explicitly. OTFP internals use C-style indices in
  the C API; Tcl scripts often translate from user-facing 1-based data.

The initial Fortran layer should be thin. Its job is to map GEMS data
structures into the C facade and return forces/energies to GEMS.

## Work To Do Before Asking Another Agent To Integrate

1. Keep `docs/developer-architecture.md` current with core/adapter boundaries.
2. Add or finalize a C embedding facade for the engine lifecycle.
3. Add a minimal C embedding example that does not source Tcl.
4. Add a minimal Fortran `iso_c_binding` example or pseudocode page.
5. Add a test that exercises the embedding path without the mock Tcl adapter.
6. Decide how OTFP will be vendored in GEMS: subtree under `lib/`, direct source
   import under a free-energy module, or a temporary development copy.
7. Only after the C/Fortran boundary is working, evaluate whether any C modules
   should be rewritten in Fortran.

## What Not To Do Yet

- Do not make GEMS call the SWIG/Tcl API.
- Do not require NAMD concepts in the engine-facing API.
- Do not translate all CV/chapeau code to Fortran before integration.
- Do not change numerical behavior while only moving adapter boundaries.

