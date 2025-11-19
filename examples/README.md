
# Examples

To run the examples, set the `CFACV_BASEDIR` environment variable to the OTFP
prefix:

    cd examples/OTFP/
    export CFACV_BASEDIR=../../usr

Examples in each folder can be run with ./run.sh. It 


# Testing the code

The `./test.sh` script run examples in test mode, which implies only a few
steps per example. It also uses the `../src/mock.tcl` script, which provides a
**mock of NAMD** for testing and development purposes. Thus, testing is
independent of NAMD distribution.

# Build examples

Build a namd script following style of other examples. The script should be
defaulted to test mode:

    set test true

and dcd output in this mode must be done in each step:

    if $test {
      dcdfreq    1
    }

Run this scirpt to build the testing files. For example:

    namd3 job0.namd 2>&1 | tee job0.log

This generate `phi.r`, `psi.r`, and similar files. Also create `xyz` file from
`dcd` with `catdcd`. For example:

    catdcd -o job0.xyz -otype xyz -s ../systems/alad.pdb job0.dcd

Then, modify `test` array in `./test.sh`  to use these files by declaring a new test.
For instance:

    [OTFP2D_periodic]="job0.namd:tol_abs:1e-4:phi.r,psi.r,job0.mk1,job0.mk2 \
          none:lambda:1e-4:job0_ch0.LAMBDA "

and run it. After fail, copy relevant files into subfolder `./reference` inside
the example path.

