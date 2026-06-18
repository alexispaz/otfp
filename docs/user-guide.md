# User guide

This page is a gentle introduction for students or new users who want to
understand what OTFP does before reading the C code.

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

## Minimal workflow

The NAMD/Tcl workflow and the GEMS embedding workflow use the same conceptual
steps:

1. Define centers.
2. Define CVs from those centers.
3. Define restraints on the CVs.
4. Optionally define chapeau/FES settings.
5. For each simulation step, copy center positions into OTFP.
6. Ask OTFP to compute CVs, bias forces, target evolution, and FES updates.
7. Copy the returned center forces back into the driver.

The exact API depends on the driver. Existing NAMD examples use Tcl scripts.
New integrations, including GEMS, should prefer the C facade documented in
[GEMS integration](integration-gems.md).

## A small mental model

It can help to think of the runtime state as five boxes:

```text
centers       static definitions of atom groups
cvs           functions of centers
restraints    bias potentials and evolving targets
center_step   positions in, forces out, for the current step
fes           chapeau grids and solve/output scheduling
```

The old code calls the whole runtime context `DataSpace`. That name is kept for
the historical C/Tcl API, but the internals are being split into the clearer
boxes shown above.

## What should a new user run first?

After building and installing into `usr/`, run the fast examples:

```sh
meson test -C build examples
```

To see the driver-facing API without NAMD, run the C embedding example:

```sh
cc examples/embedding/embedding.c \
  -Iusr/include -Lusr/lib -Wl,-rpath,"$PWD/usr/lib" -lcfacv \
  -o /tmp/otfp-embedding
/tmp/otfp-embedding
```

That example is the closest starting point for understanding how another code,
such as GEMS, should call OTFP.
