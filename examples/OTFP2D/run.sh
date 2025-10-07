#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# Runs the full example to generate meaningful results.
sed -i 's/set test true/set test false/g' job0.namd
sed -i 's/set test true/set test false/g' job1.namd

# Run simulations
namd3 job0.namd 2>&1 | tee job0.log
namd3 job1.namd 2>&1 | tee job1.log

# Go back
sed -i 's/set test false/set test true/g' job0.namd
sed -i 's/set test false/set test true/g' job1.namd

for file in $(ls *bsp); do

  f=${file%.*}

  # Decoding the FES
  ../../src/catbinsp -f $file -ol 1 \
    | awk 'BEGIN{a=""} (NR>1&&$2!=a){a=$2;$1="";$2="";print $0}'\
    > ${f}.LAMBDA

  # Unfolding the FES
  awk -v Nx=101 -v Ny=101 \
      -v xmin=2.5 -v ymin=2 \
      -v xmax=5.3 -v ymax=10 \
      'END{
        aux=(xmax-xmin)/(Nx-1)
        auy=(ymax-ymin)/(Ny-1)
        for (j=0;j<Ny;j++){
          for (i=0;i<Nx;i++){
            a=i+j*Nx+1
            print xmin+aux*i,ymin+auy*j,$a
          }
          print ""
        }
      }' ${f}.LAMBDA > ${f}.fes

  # Plotting the FES
  gnuplot << HEREGNUPLOT

  set terminal png;
  set size square
  set encoding utf8

  Is0(x)=((x==0.) ? 1/0 : x)

  stats '${f}.fes' u (Is0(\$3))
  fesmin=STATS_min

  set cbrange [0:25]
  set output "${f}_fes.png"
  set ylabel 'H' norotate
  set xlabel 'C' offset 1,0
  set pm3d map
  splot [2.5:4.5][2.5:6]\
      '${f}.fes' u 1:2:(Is0(\$3)-fesmin) notit

HEREGNUPLOT

done

#Final message
echo "Test completed sucesfully, compare the plots obtained with the reference"

