#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# Runs the full example to generate meaningful results.
sed -i 's/set test true/set test false/g' job0.namd

# Run simulations
namd3 job0.namd 2>&1 | tee job0.log

# Go back
sed -i 's/set test false/set test true/g' job0.namd

# echo "skiping input generation..."
# rm -r input/
# mkdir input/
# vmd -e -dispdev text input.tcl

for f in input/frame_*.pdb; do

    name=$(basename $f)
    name=${name%.pdb}

    sed -e "s|@INPUT|$f|" job0.namd > tmp.namd
    
		echo $f
		$namd +p2 tmp.namd > job0.log 2>&1 

		mv job0.log output/$name.log
		mv job0.coor output/$name.coor
		mv phi.r output/$name.phi.r
		mv psi.r output/$name.psi.r
		mv phi.x output/$name.phi.x
		mv psi.x output/$name.psi.x

done

# echo "skiping output.xyz generation..."
# catdcd -o reference/output.xyz -otype xyz -s input/frame_00001.pdb -namdbin output/frame_*.coor                         

#Final message
echo "Test completed sucesfully, compare the plots obtained with the reference"

