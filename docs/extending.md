# Extending OTFP

This section is for developers or integration agents who need to change OTFP,
embed it in another program, or continue the refactor.

## Current Direction

OTFP is being prepared as a vendored engine for GEMS while preserving the
existing NAMD/Tcl workflow. New core code should avoid depending on Tcl, SWIG,
or NAMD concepts. Those belong at adapter boundaries.

The preferred external C boundary is `otfp_engine.h`. Lower-level
`DataSpace_*` functions are part of the historical C/Tcl API and should not be
the default integration target for new drivers.

## What To Read

1. [Architecture](architecture.md) explains the current code layers, runtime
   objects, and ownership boundaries.
2. [GEMS integration](integration-gems.md) explains how a Fortran driver should
   call the C facade.
3. [Testing](testing.md) lists the fast guardrails that should run after C
   refactors.
4. [Code audit](audit.md) records current simplification targets and known
   risks.

## Refactor Rules

- Preserve scientific behavior unless a feature is explicitly removed and
  documented.
- Prefer mechanical splits before semantic rewrites.
- Keep fast tests passing while moving code.
- Update user and developer documentation in the same change as visible
  behavior, API, or architecture changes.
- Keep the Tcl/SWIG adapter narrow. Do not expose internal helpers through Tcl
  unless the supported workflow requires them.

## Good Extension Targets

- Add focused tests for geometry routines, chapeau restart/load paths, and
  repeated Fortran embedding steps.
- Continue splitting chapeau code by lifecycle, IO, solver, evaluation, and
  accumulation responsibilities.
- Keep moving core behavior behind clear objects: center definitions, CVs,
  restraints, step state, and FES/chapeau state.
- Replace hard exits with status returns at the C boundary, then translate
  those failures into Tcl errors in the Tcl adapter.
