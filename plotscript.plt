#!/usr/bin/gnuplot -persist
#
set key center top title " "
set ylabel "energy [a.u.]"
set y2label "dipole moment"
set grid x y
set ytics nomirror
set tics out
set y2tics
set autoscale  y
#set autoscale y2
plot 'tool.data' using 1 with linespoints t 'energy' axis x1y1, 'tool.data' using 2 t 'cell dipole moment' lt -1 axis x2y2 
#plot 'tool.data' using 1:2
pause 5
reread  
