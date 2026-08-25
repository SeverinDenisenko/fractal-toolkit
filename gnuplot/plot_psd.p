set logscale xy
set xlabel "f"
set ylabel "P(f)"
set title "PSD"
set datafile commentschars "#"
plot "psd.dat" using 1:2 with points pt 7 ps 0.5 title "P(f)"
pause -1 "Press Enter to close"
