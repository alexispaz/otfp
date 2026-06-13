# Features

This page describes the functional surface that should be preserved or
intentionally replaced during refactoring.

## NAMD/TclForce Integration

OTFP is loaded from a NAMD input script through:

```tcl
tclforcesscript $CFACV_BASEDIR/include/calcforces.tcl
source $CFACV_BASEDIR/include/cfacv.tcl
```

Runtime flow:

1. `cfacv.tcl` loads `libcfacv.so` and the Tcl helper files.
2. `parse_pdb` reads `labelPDB`; atom occupancy stores mass and beta stores
   the 1-based center id.
3. `NewDataSpace` allocates the C-side state.
4. `cvs_setup` registers collective variables.
5. `rlist_setup` registers restraints.
6. `chapeau_setup` registers OTFP/FES state when `FES(...)` is defined.
7. NAMD calls `calcforces` every step.
8. `Tcl_UpdateDataSpace` copies group positions into C, computes CVs,
   computes restraint forces, copies forces back to NAMD groups, and adds
   restraint energy.

The test runner uses `src/mock.tcl` to simulate enough NAMD behavior for fast
regression tests.

## Collective Variables

Configured with:

```tcl
set cv(num) 1
set cv(1.type) BOND
set cv(1.centers) {1 2}
```

Supported C CV types:

| Type | Centers | Notes |
| --- | ---: | --- |
| `BOND` | 2 | Distance between two centers. |
| `HALFBOND` | 2 | Distance with force only on the first center. |
| `BONDS` | even number | Sum of pair distances. |
| `ANGLE` | 3 | Angle over three centers. |
| `DIHED` | 4 | Dihedral angle. |
| `S` | even number | Bond-network switching CV. |
| `CARTESIAN_X` | 1 | X coordinate of one center. |
| `CARTESIAN_Y` | 1 | Y coordinate of one center. |
| `CARTESIAN_Z` | 1 | Z coordinate of one center. |
| `COGX` | one or more | Average x coordinate of centers. |
| `COGY` | one or more | Average y coordinate of centers. |
| `COGZ` | one or more | Average z coordinate of centers. |
| `LINE` | one or more | Projection between two reference PDB states. |
| `RMSD` | one or more | Mean squared displacement from a reference PDB. |
| `CHARMMD` | 4 | CHARMM-style dihedral energy CV using `n`, `k`, and `d`. |
| `GAUSS2D` | 5 | Two-dihedral Gaussian function. |
| `D2CHAP` | 5 | Two-dihedral CV evaluated from a chapeau file. |
| `ZSDCIRCLE` | one or more | Z-spread in a cylindrical selection. Uses global setup state. |
| `ZSDXRANGE` | one or more | Z-spread in an x range. Uses global setup state. |
| `ZSDRING` | one or more | Z-spread in a ring. Uses global setup state. |

CV output:

```tcl
set cv(1.outfile) "phi.x"
set cv(1.outfreq) 100
```

CV boundaries:

| Boundary | Behavior |
| --- | --- |
| `NADA` | No CV boundary force. |
| `SOFTUPPER` | Harmonic upper wall using `max` and `boundk`. |
| `SOFTLOWER` | Harmonic lower wall using `min` and `boundk`. |
| `SOFT` | Both soft walls using `min`, `max`, and `boundk`. |
| `AMD` | Adds the CV value as a local bias and applies unit force to the CV. |

## Restraints

Configured with:

```tcl
set restr(num) 1
set restr(1.cv) 1
set restr(1.k) 1500
set restr(1.type) tamd
```

A restraint is a linear combination of CVs. If `restr(i.cvc)` is omitted,
equal weights are assigned across `restr(i.cv)`.

Restraint potential functions:

| Function | Behavior |
| --- | --- |
| `HARMONIC` | Harmonic potential. |
| `HARMCUTO` | Harmonic potential with cutoff behavior. |

Restraint boundary functions:

| Boundary | Behavior |
| --- | --- |
| `NADA` | No restraint boundary. |
| `SOFTUPPER` | Harmonic upper wall. |
| `SOFTLOWER` | Harmonic lower wall. |
| `SOFT` | Both soft walls. |
| `PBC` | Periodic wrapping over configured `min`/`max`. |
| `PERIODIC` | Periodic wrapping over `[-pi, pi]`. |

Restraint evolution types:

| Type | Required parameters | Purpose |
| --- | --- | --- |
| `TAMD` | `g`, `dt`, plus `temp` or `biasf` | Brownian evolution of the restraint target. |
| `SMD` | `target`, `ini`, `end` | Constant-velocity steering. |
| `ADAM` | `a`, `b1`, `b2`, `e`; optional `tmin`, `decay` | Adaptive target evolution with termination. |

Restraint output:

```tcl
set restr(1.outfile) "phi.r"
set restr(1.outfreq) 100
```

## OTFP/FES Reconstruction

Configured by defining `FES(...)` before sourcing `cfacv.tcl`:

```tcl
set FES(dm) 2
set FES(min) {-3.141592653589793 -3.141592653589793}
set FES(max) { 3.141592653589793  3.141592653589793}
set FES(knots) {151 151}
set FES(period) {1 1}
set FES(start) 0
set FES(solve.start) 0
set FES(solve.freq) 2500
set FES(outfile) job0
set FES(outfreq) 2500
set FES(outlvl) 3
```

Capabilities:

- One-dimensional and two-dimensional chapeau grids.
- Optional periodic grid dimensions.
- Accumulation of force statistics from TAMD restraints.
- Dense LAPACK solve path.
- Binary output via `*_chN.bsp`.
- Restart state via `*.chN`.
- Lambda export through `catbinsp`.
- Serialization helpers for replica-exchange communication.

## Restarts

Supported restart data:

- Restraint targets through `ds_saverestrains` and `ds_loadrestrains`.
- Chapeau state through `chapeau_savestate`, `chapeau_loadstate`, and
  `chapeau_loadlambda`.
- Example restart scripts set `restart_otfp` or `restart_chapeau`, then source
  the same workflow.

## Replica Workflows

The scripts under `scripts/RE*.namd` support NAMD replica workflows:

- `REMD_NVT.namd` and `REMD_NPT.namd`: replica exchange molecular dynamics.
- `REOTFP_NVT.namd` and `REOTFP_NPT.namd`: replica exchange with OTFP.
- `REOFP_NPT.namd`: replica exchange using OTFP free-energy estimates.

Replica OTFP features include:

- `NUMREP` chapeau objects.
- Per-replica `chid` remapping on restraints.
- Chapeau serialization and aggregation across replicas.
- Restart script generation.

## Command-Line Tools

| Tool | Purpose |
| --- | --- |
| `catbinsp` | Read binary chapeau output and print selected output levels. |
| `chapsim` | Symmetrize two-dimensional chapeau state. |
| `chapconv` | Convert old chapeau state format to current format. |
| `chapadd` | Add compatible chapeau states. |
| `chapcrop` | Crop or expand a chapeau grid. |

## Fast Test Coverage

The fast regression suite currently covers mock-driven examples:

- `OTFP`
- `OTFP2D`
- `OTFP2D_periodic`
- `adam`

Run it through Meson:

```sh
meson test -C build examples
```
