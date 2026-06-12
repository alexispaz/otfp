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

The script runs from `examples/` and checks selected output files against each
example's `reference/` directory.

## Current coverage

The current runner focuses on high-level example behavior. It is valuable
because it exercises Tcl wrappers, the mock NAMD interface, chapeau output, and
numeric regression files together.

Recommended next additions:

- Small C unit tests for `measurements.c` geometry routines.
- Chapeau serialization/load round-trip tests.
- Error-path tests for missing or incompatible restart files.
- A CI job that builds, installs into `usr/`, and runs `meson test`.
