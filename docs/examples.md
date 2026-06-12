# Examples

Examples live under `examples/`. The most useful regression examples are:

| Example | Description |
| --- | --- |
| `OTFP` | One-dimensional OTFP workflow. |
| `OTFP2D` | Two-dimensional OTFP workflow. |
| `OTFP2D_periodic` | Two-dimensional periodic collective variable workflow. |
| `adam` | Adaptive-moment style workflow. |

Most runnable examples use:

```sh
cd examples/OTFP
export CFACV_BASEDIR=../../usr
./run.sh
```

The regression runner uses `src/mock.tcl`, so it can exercise Tcl-level logic
without a NAMD binary.

## Adding an example test

1. Keep the NAMD input in test mode by default.
2. Generate compact deterministic outputs.
3. Add reference files under `examples/<name>/reference/`.
4. Register the output checks in `examples/test.sh`.
5. Run `meson test -C build examples` after installing into `usr/`.

Use absolute tolerances only where numeric noise is expected; otherwise prefer
exact comparisons.
