#!/usr/bin/gnuplot -persist

# termianl
set terminal pdf enhanced color font 'Helvetica, 20' size 21 cm, 14.8 cm dashed

# picture
set grid
set xrange [ pi/200 : 1.5*pi ]
set logscale xy
set xlabel "q" 
set key box notitle right bottom spacing 1.2

set output sprintf("all.pdf")
set title "Inverse Green"
set ylabel "G^{-1}"
set yrange [ * : * ]

list=system('ls -1B eta/N=*')
plot 16*sin(x/2)**4 lw 4 lt 2 lc rgb "red", .55*x**3.3 lw 4 lt 2 lc rgb "red", for [file in list] file pt 7 ps .4 t file

# fit
p8=.30
file='eta/fit'
set output sprintf("fit.pdf")
set fit errorvariables
a=1.2
h=.75
fit [ 0 : .1 ] x**4*(a*p8/x)**h file via a, h
set xrange [ pi/150 : pi/2 ]
set label sprintf("eta\t= %.2f +/- %.2f", h, h_err) at .102,.6
set label sprintf("a\t= %.2f +/- %.2f", a, a_err) at .102,.2
set label sprintf("<-- fit region") at .098,.0007 right
set title sprintf("Inverse Green (p8 = %.2f)", p8)
plot 16*sin(x/2)**4 lw 4 lt 2 lc rgb "red", x**4*(a*p8/x)**h lw 4 lt 2 lc rgb "red", file pt 7 ps .4 t file

#    EOF
