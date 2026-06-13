# OTFP

OTFP computes free-energy surfaces in molecular dynamics simulations through
temperature-accelerated on-the-fly parameterization.

The codebase provides:

- C libraries exposed to Tcl through SWIG.
- Tcl scripts intended for NAMD TclForce workflows.
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
