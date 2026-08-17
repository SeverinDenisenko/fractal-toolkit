set logscale xy
set xlabel "f"
set ylabel "P(f)"
set title "Berg PSD"
set datafile commentschars "#"
plot "output.dat" using 1:2 with points pt 7 ps 0.5 title "P(f)"
pause -1 "Press Enter to close"
