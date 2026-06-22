# User Guide

This section is for users who need to build OTFP, install it into a local
prefix, and run the examples. It assumes basic command-line familiarity, but it
does not assume that the reader has worked with NAMD TclForces or the OTFP C
code before.

## What problem does OTFP solve?

Molecular dynamics gives atom positions and forces as a function of time. Many
questions, however, are easier to ask in terms of a few collective variables
(CVs), such as a distance, an angle, a dihedral, or a reaction coordinate. OTFP
uses those CVs to estimate a free-energy surface while the simulation is
running.

In practical terms, OTFP sits between a molecular dynamics driver and a biasing
method:

1. The driver provides the current positions of selected centers.
2. OTFP computes CV values from those centers.
3. Restraints apply bias forces on the CVs.
4. OTFP converts those CV forces back into forces on the centers.
5. When FES reconstruction is enabled, chapeau basis functions accumulate the
   statistics needed to reconstruct the free-energy surface.

## Vocabulary

| Term | Meaning in this code |
| --- | --- |
| Center | A group of one or more atoms represented by a weighted position, usually a center of mass. |
| CV | A collective variable computed from one or more centers. Examples: distance, angle, dihedral. |
| Restraint | A bias potential applied to one CV or a linear combination of CVs. |
| Target | The restraint reference value. It can be fixed or evolved by TAMD, SMD, or ADAM. |
| Chapeau | The basis/kernel function used to accumulate free-energy statistics on a grid. |
| FES | Free-energy surface reconstructed from accumulated statistics. |

## What should a new user run first?

Start by building and installing OTFP into `usr/` inside the repository:

```sh
meson setup build --prefix="$PWD/usr"
meson compile -C build
meson install -C build
```

Then run the fast example suite:

```sh
meson test -C build examples
```

The examples use a mock Tcl/NAMD adapter, so they are useful for checking the
OTFP setup quickly without launching a full NAMD calculation.

## Reading the examples

Most user-facing examples are Tcl scripts that configure centers, collective
variables, restraints, and optional FES reconstruction before the driver starts
the force loop. Read them in this order:

1. `docs/build.md` for compiler, Meson, and installation steps.
2. `docs/examples.md` for runnable examples and the mock regression workflow.
3. `docs/tcl-reference.md` for the Tcl variables used by OTFP.
4. `docs/features.md` for the complete list of supported CVs, restraints,
   chapeau/FES settings, restarts, and tools.

The NAMD/Tcl workflow follows these conceptual steps:

1. Define centers.
2. Define CVs from those centers.
3. Define restraints on the CVs.
4. Optionally define chapeau/FES settings.
5. For each simulation step, copy center positions into OTFP.
6. Ask OTFP to compute CVs, bias forces, target evolution, and FES updates.
7. Copy the returned center forces back into the driver.

To see the driver-facing API without Tcl or NAMD, run the C embedding example:

```sh
cc examples/embedding/embedding.c \
  -Iusr/include -Lusr/lib -Wl,-rpath,"$PWD/usr/lib" -lcfacv \
  -o /tmp/otfp-embedding
/tmp/otfp-embedding
```

That example is the closest starting point for understanding how another code,
such as GEMS, should call OTFP.
