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

p8=.30
a=1.3
h=.75

list=system('ls -1B eta/N=*')
plot (2*sin(x/2))**4 lw 4 lt 2 lc rgb "red", (sqrt(8)*sin(x/sqrt(8)))**4 lw 4 lt 2 lc rgb "red", x**4*(a*p8/x)**h lw 4 lt 2 lc rgb "red" t sprintf("x**4*(%.1f*p8/x)**%.2f",a,h), for [file in list] file pt 7 ps .3 t file

# fit
set fit quiet
set fit logfile '/dev/null'
file='eta/fit'
set output sprintf("fit.pdf")
set fit errorvariables
f(x) = x**4*(a*p8/x)**h
fit2 = 0.13
fit1 = .5*fit2
fit [ fit1 : fit2 ] x**4*(a*p8/x)**h file via a, h
#print "\neta=",h
#print "a  =",a,"\n"

set xrange [ pi/60 : 1 ]
set label 1 sprintf("eta\t= %.2f +/- %.2f", h, 3*h_err) at .102,.6
set label 2 sprintf("a\t= %.2f +/- %.2f", a, 3*a_err) at .102,.2
set label 3 sprintf("---|") at fit2,.0005 right
set label 4 sprintf("|---") at fit1,.0005 left
set title sprintf("Inverse Green (p8 = %.2f)", p8)
plot x**4 lw 4 lt 2 lc rgb "red", f(x) lw 4 lt 2 lc rgb "red", file  pt 7 ps .3 t file 

fit1=0.1
fit2=0.35
m=1000
step=(fit2-fit1)/m
fitcmd(i) = sprintf("fit [ .7*(fit1+step*%d) : fit1+step*%d ] x**4*(a*p8/x)**h file via a, h", i, i)
s=0
do for [i=0:m] {
    eval(fitcmd(i))
    eval("print (fit1+step*i), h")
    eval("s = s+h")
}

# plot
set terminal pdf enhanced color font 'Helvetica, 20' size 21 cm, 14.8 cm dashed
set output sprintf("fita.pdf")
set title "effective eta"
set ylabel "eta"
set xrange [ * : * ]
unset logscale xy
set key box notitle right top spacing 1.2
unset label 1
file='fita.dat'
plot file w l lw 4 t file



#    EOF
