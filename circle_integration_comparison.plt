w = 'circle_conv_data.txt'

set terminal qt
set y2tics 0, 0.1
set ytics nomirror

plot w using 1:(100*abs($2)) with lines title 'Percent Error' axis x1y1, w using 1:3 with lines 'Calculated Value' axis x1y2

set title 'Polar Integral Method'
set xlabel 'Grid size, n (nxn)'

