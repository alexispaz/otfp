# Build

OTFP can be built with Meson. A typical local build installs into `usr/`
inside the repository:

```sh
meson setup build --prefix="$PWD/usr"
meson compile -C build
meson install -C build
```

The examples expect the installed Tcl scripts and shared libraries under that
prefix.

The installed `otfp_engine.h` header is the supported C embedding facade for
external drivers. The lower-level C headers are still internal refactor
surfaces unless documented otherwise.

## Dependencies

Required build tools and libraries:

- C compiler.
- Meson and Ninja.
- SWIG.
- Tcl development package.
- BLAS, LAPACK, and LAPACKE with pkg-config metadata.

The BLAS/LAPACK provider can be selected at setup time:

```sh
meson setup build --prefix="$PWD/usr" \
  -Dblas=blas \
  -Dlapack=lapack \
  -Dlapacke=lapacke
```

## Debug builds

For debugging with NAMD or the mock runner:

```sh
meson setup build --prefix="$PWD/usr" --buildtype=debug \
  -Dc_args='-g3 -O0 -fno-omit-frame-pointer' \
  -Dc_std=gnu11
```

For sanitizer builds:

```sh
meson setup build --prefix="$PWD/usr" --buildtype=debug \
  -Db_sanitize=address,undefined \
  -Dc_args='-g3 -O0 -fno-omit-frame-pointer' \
  -Dc_std=gnu11
```
