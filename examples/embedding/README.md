# C Embedding Example

This example shows the driver-facing C API without Tcl, SWIG, NAMD, or the
mock adapter. It creates a two-center engine, registers a `BOND` CV, adds a
harmonic restraint with SMD target evolution, computes one step, and prints the
bias energy and center forces.

Build OTFP first:

```sh
meson setup build --prefix="$PWD/usr" --reconfigure
meson compile -C build
meson install -C build
```

Compile and run the example from the repository root:

```sh
cc examples/embedding/embedding.c \
  -Iusr/include -Lusr/lib -Wl,-rpath,"$PWD/usr/lib" -lcfacv \
  -o /tmp/otfp-embedding
/tmp/otfp-embedding
```

GEMS should follow this path conceptually: include `otfp_engine.h`, keep the
engine handle opaque, pass coordinates each step, then copy the returned forces
and energy into the driver.
