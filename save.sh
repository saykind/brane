#!/bin/bash

rm -r eta
rm pr.dat
for N in  100 110 120 130 140 150 160 170 180  
do
	./a.out M=0 N=$N
done
./plot.gp &> fita.dat
