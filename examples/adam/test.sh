#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

export CFACV_BASEDIR="../../usr"

namd=namd3

rm -rf output/
mkdir output/ 

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
 
#Final message
echo "Test completed sucesfully, compare the plots obtained with the reference"

