# Code Audit

This project has a mixed C/Tcl/NAMD workflow, so unused-code detection should be
treated as a candidate list rather than an automatic deletion list.

## Static candidates

Use:

```sh
scripts/audit-unused-c-symbols.sh
```

The script lists C functions with few references in `src/`. Symbols exported
through SWIG, Tcl, headers, or command line tools may still be public API and
need manual review before removal. Increase or reduce the threshold with:

```sh
MAX_REFS=3 scripts/audit-unused-c-symbols.sh
```

## Review policy

Classify each candidate as one of:

- `public-api`: exported through SWIG, Tcl, headers, or user workflows.
- `tool-entry`: used by a command-line executable or external script.
- `dead`: no internal or documented external use.
- `defer`: unclear usage, keep until covered by tests.

Only remove `dead` routines after the example regression suite passes.

## Areas to inspect first

- Old chapeau conversion compatibility code.
- Collective variable variants in `src/cvs.c`.
- Experimental or generated files under `examples/`, `scripts/ana_tamd/`, and
  `lumo/`.
