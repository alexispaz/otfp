# OTFP

OTFP computes free-energy surfaces in molecular dynamics simulations through
temperature-accelerated on-the-fly parameterization.

The codebase provides:

- A C core engine for CVs, restraints, chapeau accumulation, and bias forces.
- Tcl scripts intended for NAMD TclForce workflows.
- SWIG/Tcl adapter headers for the NAMD workflow.
- Small command-line utilities for chapeau object files.
- Mock-driven examples that can be used as regression tests without NAMD.

## Main components

| Component | Purpose |
| --- | --- |
| `libcfacv` | Collective variables, restraints, chapeau accumulation, and Tcl interface. |
| `libcvs` | Collective variable calculations and Tcl interface. |
| `catbinsp` | Export chapeau binary state data. |
| `chapsim` | Symmetrize chapeau data. |
| `chapadd` | Add compatible chapeau files. |
| `chapcrop` | Crop or expand a chapeau grid. |

The original background and citations are kept in the repository
`README.md`; this documentation focuses on build, testing, examples, and code
maintenance.

## Integration direction

OTFP is being prepared so its C core can be used as an engine by other
molecular drivers. The immediate target is integration into GEMS as vendored
source driven from Fortran/C, while preserving the NAMD/Tcl adapter as a
regression and compatibility path. See [GEMS integration](integration-gems.md)
for the current plan.
