reset

set linetype 1 lc rgb '#1f77b4'
set linetype 2 lc rgb '#ff7f0e'
set linetype 3 lc rgb '#2ca02c'
set linetype 4 lc rgb '#d62728'
set linetype 5 lc rgb '#9467bd'
set linetype 6 lc rgb '#8c564b'
set linetype 7 lc rgb '#e377c2'
set linetype 8 lc rgb '#7f7f7f'
set linetype 9 lc rgb '#bcbd22'
set linetype 10 lc rgb '#17becf'
set linetype cycle 10

set term pngcairo size 1000, 1200 font "Open Sans, 20"

set output "vertical_sensitivity.png"
set xlabel "Relative Jacobian q_{SO_2}(dJ/dq_{SO_2})"
set ylabel "Altitude (km)"
set grid
set ytics 5
set mytics 5
set key bottom
plot "vertical_sensitivity.dat" u 4:($1*1e-3) w l lw 3 title "SO_2"