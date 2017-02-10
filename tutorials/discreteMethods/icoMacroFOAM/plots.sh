gnuplot -p -e 'set xlabel "Iteration"; set ylabel "mass flow rate (10^{-12} Kg/s)";set yrange[0:6]; p "multiscale/mDots_vs_iteration.xy" u 1:($2*1e12) w l, "" u 1:($3*1e12) w l, "" u 1:($4*1e12) w l, "" u 1:($5*1e12) w l, "" u 1:($6*1e12) w l'

#gnuplot -p -e 'set xlabel "Iteration"; set ylabel "force (10^{-13} N)";set yrange[0:10]; p "multiscale/forces_vs_s.xy" u 1:($2*1e13) w l, "" u 1:($3*1e13) w l, "" u 1:($4*1e13) w l, "" u 1:($5*1e13) w l, "" u 1:($6*1e13) w l'
