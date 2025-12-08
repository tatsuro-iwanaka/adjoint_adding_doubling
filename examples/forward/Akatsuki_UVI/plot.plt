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

set term pngcairo size 1200, 800 font "Open Sans, 20"

set output "reflectrance.png"
set xlabel "Emission zenith angle (degree)"
set ylabel "Reflectance"
set xrange[0:90]
set grid
plot "venus_35.dat" w lp lw 3 pt 7 ps 0.5 title "d=35", "venus_50.dat" w lp lw 3 pt 7 ps 0.5 title "d=50"


set output "thermal_emission.png"
set xlabel "Emission zenith angle (degree)"
set ylabel "Radiance (W/m^2/sr/m)"
set xrange[0:90]
set grid
plot "venus_35.dat" u 1:3 w lp lw 3 pt 7 ps 0.5 title "d=35", "venus_50.dat" u 1:3 w lp lw 3 pt 7 ps 0.5 title "d=50"
