# Benchmarks

Cost of OTFP in a production-size NAMD run, measured on the CCAD (UNC)
cluster *mendieta* between 2026-10-04 and 2026-10-08.

## Setup

- **System:** SOD1, 217,291 atoms, PME, 2 fs time step, harmonic CA restraints.
- **NAMD:** 3.0.3 built with Tcl (`interface=tcl`), so TclForces is available.
  A build without Tcl was used as a baseline for the Tcl build itself.
- **OTFP:** `master` at commit `8019bb3`, loaded through TclForces as described
  in the [feature reference](features.md).
- **Hardware:** nodes with 2x Intel E5-2680 v2 (20 cores) and 2x NVIDIA A30,
  ConnectX-3 QDR InfiniBand.
- **Modes:**
    - *GPU-resident* (`+single_node_gpu`, charm++ multicore), 1 GPU.
    - *GPU offload* over UCX (charm++ `ucx`, PMIx, SMP), one process per GPU,
      1 to 4 nodes.
- **Metric:** ns/day reported by NAMD. Numbers in parentheses are Slurm job
  ids on mendieta.

## Results

### GPU-resident, 1 A30

| run | OTFP off | OTFP on | OTFP cost |
|---|---|---|---|
| repeat, run 1 | 27.66 (1179695) | 19.69 (1179696) | |
| repeat, run 2 | 27.68 (1179702) | 19.72 (1179703) | |
| first round | 27.77 (1179398) | 19.60 (1179400) | |
| **summary** | **27.7** | **19.7** | **-29 %** |

The spread between runs is below 1 %.

The Tcl build by itself costs nothing: without OTFP, the Tcl build gave
27.77 and 27.81 ns/day (1179398, 1179410) and the build without Tcl 27.93 and
27.69 (1179402, 1179414), on two software stacks built with different
compilers.

### GPU offload (UCX), first round, one run per point

| stack | nodes | OTFP off | OTFP on |
|---|---|---|---|
| A | 1 (2 GPUs) | 15.00 (1179404) | 14.93 (1179407) |
| A | 2 | 17.39 (1179506) | 17.08 (1179508) |
| A | 4 | 17.53 (1179554) | 19.20 (1179509) |
| B | 1 (2 GPUs) | 15.04 (1179416) | 14.83 (1179419) |
| B | 2 | 17.24 (1179510) | 15.96 (1179512) |
| B | 4 | 19.56 (1179511) | 18.84 (1179513) |

Stacks A and B are the same NAMD 3.0.3 + OTFP built with GCC 12.5 and 11.2.
Differences of up to about 10 % between single runs are noise on this
cluster (the 4-node pair of stack A is even inverted), so the OTFP cost in
offload mode, 0-8 %, is not distinguishable from zero.

Older builds in offload mode gave the same picture: NAMD 2.14 with OTFP
`a0c905e` and NAMD 3.0b6 with OTFP `adam` were 0.1-1.1 ns/day slower with
OTFP on 1 and 2 nodes (1179422-1179425, 1179514-1179517).

## Reading

- OTFP costs about 29 % in GPU-resident mode and is within noise in offload
  mode.
- NAMD calls `calcforces` every step, and NAMD 3.0.3 has no keyword to run
  TclForces/GlobalMaster every N steps. A likely cause of the GPU-resident
  cost is that each call brings the center coordinates from the GPU to the
  host; in offload mode they are already there. This has not been profiled.
- Even with OTFP, GPU-resident on one A30 (19.7 ns/day) is the best per-GPU
  option, about equal to offload on 4 nodes (18.8-19.2 ns/day).
- GPU-resident on 2 GPUs cannot run this input: NAMD stops with "Restraints
  (constraints on) are not compatible with multi-GPU GPUresident".
