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
    cleanup.
  - `cfacv_restraint.c`: restraint construction, restart state, energy,
    boundaries, and target evolution.
  - `cfacv.c`: force-loop orchestration and reporting.
  - Future chapeau splits should separate IO/serialization, solver, evaluation,
    and accumulation.
- Prefer descriptor tables over large factory switches for extensible
  registries such as CV types.
- Add destructors and ownership tests before moving code that allocates memory.

## Build And Test Guardrails

Run these after each meaningful C refactor:

```sh
meson setup build --prefix=$PWD/usr --reconfigure
meson compile -C build
meson test -C build lifecycle
meson install -C build
meson test -C build examples
```

Run documentation validation after docs or architecture changes:

```sh
mkdocs build --strict --site-dir /tmp/otfp-mkdocs-site
```

The `examples` test depends on the installed local prefix under `usr/`, so run
`meson install -C build` before it.

## Documentation

- Keep user-facing features documented in `docs/features.md`.
- Keep developer architecture and refactor status in
  `docs/developer-architecture.md`.
- Keep GEMS and driver-embedding context in `docs/integration-gems.md`.
- Keep test commands and coverage notes in `docs/testing.md`.
- If a user-visible configuration key, feature, tool, or workflow changes,
  update docs in the same commit.

## Current Fast Tests

- `meson test -C build lifecycle`: small C constructor/destructor guardrail.
- `meson test -C build embedding`: C driver-embedding smoke test that avoids
  Tcl and the NAMD mock adapter.
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
  code depend on Tcl, SWIG, or NAMD concepts.
- Sparse/MKL chapeau support exists in source but is marked under construction
  in the build logic.
- Restart files serialize raw structs in places; be careful when changing
  struct layout or ownership fields.
