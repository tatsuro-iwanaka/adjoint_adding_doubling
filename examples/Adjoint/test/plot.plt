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

set output "sensitivity.png"
set title "Sensitivities on Optical Thickness and Single Scattering Albedo"
set xlabel "dJ/dx"
set ylabel "layer"
set grid
set xrange[-0.02:0.05]
set key bottom
plot "layer_sensitivity.dat" u 2:1 w l lw 3 title "optical thickness", "layer_sensitivity.dat" u 3:1 w l lw 3 title "single scattering albedo"

set output "sensitivity_phase_function.png"
set title "Sensitivity on Scattering Phase Function"
set xlabel "dJ/dx"
set ylabel "layer"
set grid
unset xrange
set key bottom
set logscale x
plot "layer_sensitivity.dat" u 4:1 w l lw 3 title "{/Symbol Q}=0", "layer_sensitivity.dat" u 5:1 w l lw 3 title "{/Symbol Q}=30", "layer_sensitivity.dat" u 6:1 w l lw 3 title "{/Symbol Q}=60", "layer_sensitivity.dat" u 7:1 w l lw 3 title "{/Symbol Q}=90", "layer_sensitivity.dat" u 8:1 w l lw 3 title "{/Symbol Q}=120", "layer_sensitivity.dat" u 9:1 w l lw 3 title "{/Symbol Q}=150", "layer_sensitivity.dat" u 10:1 w l lw 3 title "{/Symbol Q}=180"