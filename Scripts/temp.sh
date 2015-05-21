#!/bin/bash

fsta=station.lst
stainfoV=()
ista=0
while read sta lon lat; do
	stainfoV[ista]=`echo $sta $lon $lat`
	let ista++
done < $fsta
nsta=$ista

for ((ista1=0;ista1<$nsta;ista1++)); do
	lon1=`echo ${stainfoV[ista1]} | awk '{print $2}'`
	lat1=`echo ${stainfoV[ista1]} | awk '{print $3}'`
	for ((ista2=$ista1+1;ista2<$nsta;ista2++)); do
		lon2=`echo ${stainfoV[ista2]} | awk '{print $2}'`
		lat2=`echo ${stainfoV[ista2]} | awk '{print $3}'`
		dist=`/home/tianye/usr/bin/get_dist $lat1 $lon1 $lat2 $lon2 d`
		if [ `echo $dist | awk '{if($1<1){print 1}else{print 0}}'` == 0 ]; then continue; fi
		echo $dist ${stainfoV[ista1]} ${stainfoV[ista2]}
	done
done
