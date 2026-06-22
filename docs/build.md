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

Use `meson configure build` after setup to inspect available project options.
Use `meson dist -C build` when you need a source distribution.

## Documentation site

The documentation is built with MkDocs from the repository root:

```sh
mkdocs build
```

This writes the HTML site to `site/`. The directory is generated output and is
ignored by git.

For local reading while editing:

```sh
mkdocs serve
```

Then open the URL printed by MkDocs, usually `http://127.0.0.1:8000/`.

For stricter checks before committing documentation changes:

```sh
mkdocs build --strict --site-dir /tmp/otfp-mkdocs-site
```

If MkDocs fails before building pages, check that `mkdocs` and the
`readthedocs` theme are available in the active Python environment. The theme
is part of a normal MkDocs installation, but minimal or isolated environments
may need MkDocs installed explicitly.

## Dependencies

Required build tools and libraries:

- C compiler.
- Meson and Ninja.
- MkDocs for documentation builds.
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

When running NAMD with AddressSanitizer, the sanitizer runtime may need to be
preloaded:

```sh
LD_PRELOAD="$(gcc -print-file-name=libasan.so)" namd3 +p1 job0.namd
```

For optimized local builds:

```sh
meson setup build --prefix="$PWD/usr" --buildtype=release \
  -Doptimization=3 \
  -Dwarning_level=0 -Dwerror=false \
  -Db_ndebug=true -Dstrip=true
```
