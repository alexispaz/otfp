# Examples

Examples live under `examples/`. The most useful regression examples are:

| Example | Description |
| --- | --- |
| `OTFP` | One-dimensional OTFP workflow. |
| `OTFP2D` | Two-dimensional OTFP workflow. |
| `OTFP2D_periodic` | Two-dimensional periodic collective variable workflow. |
| `adam` | Adaptive-moment style workflow. |
| `embedding` | Minimal C driver example using `otfp_engine.h` without Tcl/NAMD. |

Most runnable examples use:

```sh
cd examples/OTFP
export CFACV_BASEDIR=../../usr
./run.sh
```

The regression runner uses `src/mock.tcl`, so it can exercise Tcl-level logic
without a NAMD binary.

## C embedding example

The `examples/embedding` directory demonstrates the driver-facing API intended
for GEMS and other non-NAMD programs. It includes only `otfp_engine.h`, links
against `libcfacv`, computes one step, and prints the bias energy and center
forces.

After installing locally:

```sh
cc examples/embedding/embedding.c \
  -Iusr/include -Lusr/lib -Wl,-rpath,"$PWD/usr/lib" -lcfacv \
  -o /tmp/otfp-embedding
/tmp/otfp-embedding
```

## Adding an example test

1. Keep the NAMD input in test mode by default.
2. Generate compact deterministic outputs.
3. Add reference files under `examples/<name>/reference/`.
4. Register the output checks in `examples/test.sh`.
5. Run `meson test -C build examples` after installing into `usr/`.

Use absolute tolerances only where numeric noise is expected; otherwise prefer
exact comparisons.
