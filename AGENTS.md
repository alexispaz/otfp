# Agent Instructions

This file gives repo-specific instructions for automated coding agents working
on OTFP.

## Working Branch

- Work on branch `refactor`.
- Do not make refactor commits on `master`.
- The `refactor` branch is intentionally allowed to break backward
  compatibility when that simplifies the architecture.
- Preserve scientific behavior unless a feature is explicitly removed and
  documented.

## Git Hygiene

- Commit every coherent change set.
- Use a short commit title plus a clear body explaining:
  - what changed,
  - why it changed,
  - what behavior or architecture boundary it protects.
- Do not commit generated outputs from local examples, `usr/`, `build/`, or
  ad-hoc analysis folders unless the task explicitly asks for them.
- The repository may contain many untracked local artifacts. Ignore unrelated
  untracked files.
- Never revert user or pre-existing work unless explicitly requested.

## Refactor Policy

- Prefer mechanical splits before semantic rewrites.
- Keep public behavior covered by fast tests while moving code.
- Do not preserve compatibility solely for old API shape. Remove or rename
  interfaces when the new architecture is clearer, but update docs and tests in
  the same change.
- Split large C files into ownership-focused modules:
  - `cfacv_dataspace.c`: `DataSpace` allocation, registration, accessors, and
    cleanup. `DataSpace` is a legacy runtime context, not the final conceptual
    model object.
  - `cfacv_force.c`: per-step force-loop orchestration, restraint values,
    restraint forces, CV boundary forces, FES/chapeau updates, restart
    scheduling, and target evolution.
  - `cfacv_restraint.c`: restraint construction, energy functions, boundary
    functions, and TAMD/SMD/ADAM target evolution.
  - `cfacv_restart.c`: restraint restart save/load persistence.
  - `cfacv.c`: banner, file helpers, CV computation, energy reads, and
    reporting.
  - Future chapeau splits should separate IO/serialization, solver, evaluation,
    and accumulation.
- Keep internal state names tied to the scientific/domain concept. Current
  `DataSpace` subobjects are:
  - `CenterSet`: static atom-center definitions.
  - `CvSet`: collective variables.
  - `RestraintSet`: restraints plus restraint restart policy.
  - `CenterStepState`: per-step center positions and accumulated forces.
  - `FesState`: chapeau/FES state and accumulate/solve scheduling.
- Do not reintroduce generic names such as `Registry` or ambiguous names such
  as `AtomSet` for center definitions. A center is a group/weighted position,
  not a raw atom.
- Keep `otfp_engine_*` prefixed because it is public embedding API intended for
  vendoring in GEMS. Internal `static` helpers do not need an OTFP prefix.
- Prefer descriptor tables over large factory switches for extensible
  registries such as CV types.
- Add destructors and ownership tests before moving code that allocates memory.

## Build And Test Guardrails

Run these after each meaningful C refactor:

```sh
meson setup build --prefix=$PWD/usr --reconfigure
meson compile -C build
meson test -C build lifecycle
meson test -C build embedding
meson test -C build embedding-fortran
meson test -C build restart
meson install -C build
meson test -C build examples
```

Run documentation validation after docs or architecture changes:

```sh
mkdocs build --strict --site-dir /tmp/otfp-mkdocs-site
```

The `examples` test depends on the installed local prefix under `usr/`, so run
`meson install -C build` before it.

If SWIG-visible headers change (`cfacv_tcl.h`, `cvs_tcl.h`, or a header they
include), stale generated wrappers can survive an incremental build. If the
wrapper references a removed symbol, clean and regenerate the build directory:

```sh
rm -rf build
meson setup build --prefix=$PWD/usr --reconfigure
meson compile -C build
```

## Documentation

- Keep beginner-facing explanations in `docs/user-guide.md`.
- Keep user-facing features documented in `docs/features.md`.
- Keep developer architecture and refactor status in
  `docs/developer-architecture.md`.
- Keep GEMS and driver-embedding context in `docs/integration-gems.md`.
- Keep the standalone C embedding example in `examples/embedding`.
- Keep test commands and coverage notes in `docs/testing.md`.
- Keep MkDocs usage in `docs/build.md`; generated `site/` output should remain
  ignored.
- If a user-visible configuration key, feature, tool, or workflow changes,
  update docs in the same commit.

## Current Fast Tests

- `meson test -C build lifecycle`: small C constructor/destructor guardrail.
- `meson test -C build embedding`: C driver-embedding smoke test that avoids
  Tcl and the NAMD mock adapter.
- `meson test -C build embedding-fortran`: optional Fortran `iso_c_binding`
  smoke test when a Fortran compiler is available.
- `meson test -C build restart`: C restart round-trip smoke test for restraint
  target persistence.
- `meson test -C build examples`: mock-driven Tcl/NAMD example regression
  suite.

## Code Style

- Follow the existing C style for now: compact functions, simple structs, and
  explicit allocation/free pairs.
- Use `calloc` for runtime structs that contain optional pointers.
- Close owned file handles in destructors when ownership is clear.
- Avoid broad formatting churn during mechanical moves.
- Use `rg` for code search.

## Known Constraints

- SWIG uses Tcl-specific headers (`cfacv_tcl.h`, `cvs_tcl.h`). Keep that
  adapter boundary narrow; do not expose internal helpers through Tcl unless the
  supported workflow requires them.
- OTFP is being prepared as a vendored engine for GEMS. Do not make new core
  code depend on Tcl, SWIG, or NAMD concepts. New driver-facing code should use
  `otfp_engine.h`.
- `otfp_engine.h` is the supported external C facade. Lower-level
  `DataSpace_*` functions remain the historical C/Tcl API and should not be the
  preferred integration target for GEMS.
- Sparse/MKL chapeau support exists in source but is marked under construction
  in the build logic.
- Restart files serialize raw structs in places; be careful when changing
  struct layout or ownership fields.
