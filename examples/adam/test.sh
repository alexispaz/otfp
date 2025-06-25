#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

export CFACV_BASEDIR="../../usr"

namd=namd2

# # Run simulations
$namd job0.namd 2>&1 | tee job0.log

#Final message
echo "Test completed sucesfully, compare the plots obtained with the reference"

