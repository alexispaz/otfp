set terminal png;
set size square
set encoding utf8

Is0(x)=((x==0.) ? 1/0 : x)

stats '../../OTFP2D_periodic/job0_ch0.fes' u (Is0($3))
fesmin=STATS_min

set cbrange [0:25]
set output "plot.png"
set xlabel 'φ' norotate
set ylabel 'ψ' offset 1,0
set pm3d map

# Convertir índice ($0) a color RGB (azul -> rojo)
# color(gradient) = sprintf("#%02X%02X%02X", 255*gradient, 0, 255*(1-gradient))

splot [-3.1415:3.1415][-3.1415:3.1415]\
    '../../OTFP2D_periodic/job0_ch0.fes' u 2:1:(Is0($3)-fesmin) notit w pm3d,\
    '< paste phi.r psi.r' u 1:4:(26) w p ps 0.3

