#!/bin/bash

# Adapted form GEMS: Extensible Molecular Simulator
# https://github.com/alexispaz/GEMS
# Copyright (C) 2020 Sergio Alexis Paz

set -o nounset
set -o pipefail

# Environment
export CFACV_BASEDIR=../../usr
mock=$CFACV_BASEDIR/include/mock.tcl
catb=$CFACV_BASEDIR/bin/catbinsp

# Log into named pipe
rm -rf test.log
exec > >(tee -i test.log)
exec 2>&1
 
# Timing
TIMEFORMAT="%2Rsec %2Uusr %2Ssys (%P%% cpu)."
start=$(date +%s.%N)
SECONDS=0
 
# Labels
pass='PASS'
fail='FAIL'
  
usage() {
cat <<'WXYZ'
Use:
    test.sh [directory]

Description:
  Execute `directory` test. If directory is not given, execute all tests.

Flags:
  -h this help
  -i interactive mode
  -c highlight with colors
WXYZ
exit
}

exe() {
  # valgrind --leak-check=yes --track-origins=yes ../../src/gems $1 &>> test.log
  local name="$1"; shift

  line="TEST... $name "
  echo -ne "$line \r"

  timing=$( { time tclsh $mock "$name" > /dev/null 2>&1; } 2>&1 )
	if [ $? -ne 0 ]; then
    line="ERROR ON EXECUTION"
    echo -ne "$line \r"
    return 1
  fi
	
  line="XXXX : ${timing%.} : $name \r"
  echo -ne "$line \r"
}

check() {
  local mode="$1"; shift
  local tol="$1"; shift
  local files=("$@")
  local ec=0

  for f in "${files[@]}"; do
    case $mode in

      exact)
        diff -q "$f" "reference/$f" || ec=1 ;;

      tol_abs)
        awk -v tol="$tol" '
          function abs(x){return x<0?-x:x}
          (NR==FNR){for(i=1;i<=NF;i++)a[i,NR]=$i;next}
          {for(i=1;i<=NF;i++) if(abs($i-a[i,FNR])>tol) exit 1}
        ' "$f" "reference/$f" || ec=1 ;;
				
      tol_trunc)
        awk -v fmt="$tol" '
          function trunc(x){return sprintf(fmt, x)}
          (NR==FNR){for(i=1;i<=NF;i++)a[i,NR]=$i;next}
          {for(i=1;i<=NF;i++) if(trunc($i)!=trunc(a[i,FNR])) exit 1}
        ' "$f" "reference/$f" || ec=1 ;;
    esac
  done

  [ $ec -eq 0 ] && echo -e "$line$pass" || echo -e "$line$fail"
}

run_test() {
  local name=( $(echo "$1" | tr ',' ' ') )
  local mode="$2"
  local tol="$3"
  local files=( $(echo "$4" | tr ',' ' ') ) 

  # Ejecuto
  for n in "${name[@]}"; do
		 [[ "$n" == "none" ]] && continue
     rm -f "${files[@]}"
     exe "$n"
  done 

  if [[ "$mode" == "lambda" ]]; then
	  $catb -f job0_ch0.bsp -ol 3 2>/dev/null > "${files[0]}"
		check "tol_abs" "$tol" "${files[0]}"
  else
		check "$mode" "$tol" "${files[@]}"
  fi
}

# Flags
while getopts ":hic" option; do
case $option in
   h) usage;;
   c) pass='\033[0;32mPASS\033[0m'; fail='\033[0;31mFAIL\033[0m';;
   i) prompt(){ read -p "Test $1? [Y/N] (Y) " ans; [[ $ans =~ ^[nN] ]] && return 1 || return 0; };;
   :) echo "Error: -$OPTARG requires argument"; usage;;
   ?) echo "Error: unknown option -$OPTARG"; usage;;
  esac
done
shift $((OPTIND-1))

# Prompt
if [[ -z ${1:+x} ]]; then
  prompt(){ echo ""; echo "Testing $1"; }
else
  input=${1%/}
  prompt(){ [[ $1 == $input ]]; }
fi

# Test declarations
# Format: [folder]="input:check mode:check parameter:check file"
declare -A tests=(
  [OTFP]="job0.namd:tol_abs:1e-4:job0.cv,job0.mk1,job0.mk2,job0.r \
          none:lambda:1e-4:job0_ch0.LAMBDA \
          job1.namd:tol_abs:1e-4:job1.cv,job1.mk1,job1.mk2 \
          none:tol_abs:1e-4:job1.r \
          none:lambda:1e-4:job1_ch0.LAMBDA "
  [OTFP2D]="job0.namd:tol_abs:1e-4:job0.cv1,job0.cv2,job0.mk1,job0.mk2,job0.r1,job0.r2 \
          none:lambda:1e-4:job0_ch0.LAMBDA \
          job1.namd:tol_abs:1e-4:job1.cv1,job1.cv2,job1.mk1,job1.mk2,job1.r1,job1.r2 \
          none:lambda:1e-4:job1_ch0.LAMBDA "
  [OTFP2D_periodic]="job0.namd:tol_abs:1e-4:phi.r,psi.r,phi.x,psi.x,job0.mk1,job0.mk2 \
          none:lambda:1e-4:job0_ch0.LAMBDA "
  [adam]="job0.namd:tol_abs:1e-4:phi.r,psi.r,phi.x,psi.x,job0.mk1,job0.mk2"
)
# [CVs]="rmsd,line:none:0:"

# Test execution
for section in "${!tests[@]}"; do
  if prompt "$section"; then
    cd "$section"
    # separar las entradas por espacios
    for entry in ${tests[$section]}; do
      IFS=":" read -r name mode tol files <<<"$entry"
      # expandir comas en espacios para pasar a run_test
      run_test "$name" "$mode" "$tol" "$files"
    done
    cd ..
  fi
done


# Print total time spent in the test. 
# Despite the Nanosecond format, the precision will be
# around the millisecond or probably less.
end=$(date +%s.%N)
if command -v bc &>/dev/null; then
  echo "Total time: $(echo "$end - $start" | bc -l)"
elif command -v perl &>/dev/null; then
  echo "Total time: $(perl -E "say $end-$start")"
fi

