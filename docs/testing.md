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

The script runs from `examples/` and checks selected output files against each
example's `reference/` directory.

## Current coverage

The current test set combines a low-level lifecycle test with high-level
example behavior. The example runner is valuable because it exercises Tcl
wrappers, the mock NAMD interface, chapeau output, and numeric regression files
together.

Recommended next additions:

- Small C unit tests for `measurements.c` geometry routines.
- Chapeau serialization/load round-trip tests.
- Error-path tests for missing or incompatible restart files.
- A CI job that builds, installs into `usr/`, and runs `meson test`.
