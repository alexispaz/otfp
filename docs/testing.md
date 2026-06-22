# Testing

The main regression entry point is:

```sh
examples/test.sh
```

After configuring and installing locally, it is also available through Meson:

```sh
meson setup build --prefix="$PWD/usr"
meson compile -C build
meson install -C build
meson test -C build examples
```

There is also a small C lifecycle test:

```sh
meson test -C build lifecycle
```

It constructs and frees a minimal `DataSpace` with atom centers, one CV, and
one restraint. This test is intentionally narrow; it catches constructor,
destructor, and link-level regressions before the slower example suite runs.

The C embedding facade has its own smoke test:

```sh
meson test -C build embedding
```

It constructs an `otfp_engine`, registers two centers, one `BOND` CV, and one
harmonic restraint with SMD target evolution, then computes one step without
sourcing Tcl or using the mock NAMD adapter. This is the first guardrail for
future GEMS integration.

When a Fortran compiler is available, Meson also builds and runs an
`iso_c_binding` smoke test:

```sh
meson test -C build embedding-fortran
```

This test exercises the same `otfp_engine` path from Fortran and verifies that
opaque handles, C strings, contiguous arrays, and returned force pointers can be
used from a Fortran caller.

The script runs from `examples/` and checks selected output files against each
example's `reference/` directory.

## Current coverage

The current test set combines low-level lifecycle, C embedding, and optional
Fortran embedding tests with high-level example behavior. The example runner is
valuable because it exercises Tcl wrappers, the mock NAMD interface, chapeau
output, and numeric regression files together.

Recommended next additions:

- Small C unit tests for `measurements.c` geometry routines.
- Chapeau serialization/load round-trip tests.
- Error-path tests for missing or incompatible restart files.
- More Fortran coverage for GEMS-like setup and repeated steps.
- A CI job that builds, installs into `usr/`, and runs `meson test`.
