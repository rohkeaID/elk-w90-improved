#!/bin/bash

fact=$(echo $(tail -n2 BAND.OUT |head -n1) $(tail -n2 ELK_band.dat |head -n1) | awk '{print $3/$1}')
sed "s/set nokey/#set nokey/" ELK_band.gnu | sed "s/plot \"ELK_band.dat\"/plot \"ELK_band.dat\" w l lt 2 lc rgb \"red\" title \"ELK+Wannier90\"/" 
echo " repl 'BAND.OUT' u (\$1*$fact):(\$2*27.211) w l lt rgb \"black\" title \"ELK task 20\""
echo "set ylabel 'Energy (eV)'"
echo "set term png"
echo "set output 'band.png'"
echo "repl"
