#!/usr/bin/gnuplot -persist

# termianl
set terminal pdf enhanced color font 'Helvetica, 20' size 21 cm, 14.8 cm dashed

# picture
set grid
set xrange [ pi/50 : pi ]
set logscale xy
set xlabel "q" 
set key box notitle right bottom spacing 1.2

set output sprintf("all.pdf")
set title "Inverse Green (p8=1.5)"
set ylabel "G^{-1}"
set yrange [ * : * ]


list=system('ls -1B green/N=*')
plot x**4 lw 4 lt 2 lc rgb "red", 1.5*x**3.2 lw 4 lt 2 lc rgb "red", for [file in list] file pt 7 ps .5 t file

#    EOF
