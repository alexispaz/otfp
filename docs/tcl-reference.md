# Tcl Configuration Reference

OTFP's original user interface is a Tcl layer loaded from a NAMD input script.
The Tcl files parse user variables, create the C-side objects through SWIG, and
let `calcforces.tcl` exchange positions, forces, and energies with NAMD
TclForces.

The usual entry point is:

```tcl
tclforcesscript $CFACV_BASEDIR/include/calcforces.tcl
source $CFACV_BASEDIR/include/cfacv.tcl
```

`CFACV_BASEDIR` should point to the local install prefix, for example
`../../usr` from many example directories.

## Basic Tcl Shape

The Tcl setup uses global associative arrays. Keys follow the pattern
`name(index.field)`.

```tcl
set cv(num) 1
set cv(1.type) BOND
set cv(1.centers) {1 2}

set restr(num) 1
set restr(1.cv) 1
set restr(1.k) 1500
set restr(1.type) tamd
```

This means: create one collective variable, define it as the distance between
centers 1 and 2, then create one restraint acting on that CV.

## Centers

Centers are weighted atom groups. In the NAMD workflow they are usually read
from `labelPDB`: occupancy stores the atom mass and beta stores the 1-based
center id. OTFP then registers those groups with NAMD and uses their center
positions during the force callback.

## Collective Variables

Collective variables are configured through `cv(...)`.

```tcl
set cv(num) 1
set cv(1.type) BOND
set cv(1.centers) {1 2}
set cv(1.outfile) "distance.cv"
set cv(1.outfreq) 100
```

Common types include `BOND`, `ANGLE`, `DIHED`, Cartesian coordinates, RMSD,
line projections, and chapeau-based CVs. The complete type list is in
[Feature reference](features.md).

## Restraints

Restraints are configured through `restr(...)`.

```tcl
set restr(num) 1
set restr(1.cv) 1
set restr(1.k) 1500
set restr(1.type) tamd
set restr(1.outfile) "distance.r"
set restr(1.outfreq) 100
```

A restraint can act on one CV or a linear combination of CVs. If `restr(i.cvc)`
is omitted, OTFP assigns equal weights to the listed CVs.

The main evolution modes are:

| Type | Use |
| --- | --- |
| `TAMD` | Temperature-accelerated evolution of the restraint target. |
| `SMD` | Constant-velocity steering between initial and final targets. |
| `ADAM` | Adaptive target evolution. |

## FES and Chapeau Settings

FES reconstruction is enabled by defining `FES(...)` before sourcing
`cfacv.tcl`.

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

The chapeau basis is the kernel used to accumulate force statistics on the FES
grid. Changing the basis would affect which grid points receive accumulated
statistics and may also affect the solve step.

## Restarts

Example restart scripts usually set `restart_otfp` or `restart_chapeau` before
sourcing the normal setup files. Restraint targets and chapeau state are loaded
through the Tcl/SWIG layer, then the force callback continues from the restored
state.

## Tcl References

Useful references when editing OTFP Tcl scripts:

- Tcl language manual: https://www.tcl.tk/man/tcl8.6/TclCmd/contents.htm
- NAMD TclForces documentation: https://www.ks.uiuc.edu/Research/namd/2.14/ug/node48.html
- SWIG Tcl documentation: https://www.swig.org/Doc4.1/Tcl.html
