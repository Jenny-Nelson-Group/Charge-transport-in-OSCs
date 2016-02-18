#!/usr/bin/gnuplot
reset
set terminal postscript
set output "rdf.ps"
set nokey
set ylabel "rdf"
set xlabel "r"
plot "rdf.xvg" using 1:2 with lines