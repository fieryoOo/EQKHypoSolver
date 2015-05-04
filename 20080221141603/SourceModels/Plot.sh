#!/bin/bash

types=('T' 'CORR' 'AMP')

iter=0
while [ $iter -lt 3 ]; do
 
tp=${types[$iter]}
col=2
files=(`ls *_${tp}_AZ_s206_d36_r-95_h7`)
fout=5mod_${tp}.ps

gnuplot << END
set term postscript enhanced color
set out '$fout'
plot '${files[0]}' using 1:${col} lt 1 lc rgb "red" w l,\
'${files[1]}' using 1:${col} lt 1 lc rgb "orange" w l,\
'${files[2]}' using 1:${col} lt 1 lc rgb "green" w l,\
'${files[3]}' using 1:${col} lt 1 lc rgb "blue" w l,\
'${files[4]}' using 1:${col} lt 2 lc rgb "black" w l,\
'${files[5]}' using 1:${col} lt 3 lc rgb "purple" w l
quit
END
echo $fout

let iter++
done
