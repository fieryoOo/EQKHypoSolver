#!/bin/bash

### main ###
if [ $# != 4 ] && [ $# != 5 ]; then
   echo "Usage: "$0" [type(R/L)] [iter#] [minamp] [maxamp] [amp factor (optional)]"
   exit
fi

type=$1
ampmin=$3
ampmax=$4
if [ `echo $ampmin $ampmax | awk '{if($1<0||$1>=$2){print 0}else{print 1}}'` == 0 ]; then
	echo "Invalid min||max amp: "$ampmin" "$ampmax
	exit
fi
psout=results_${type}.ps 
color=(orange lightred forestgreen steelblue midnightblue)
gmtset HEADER_FONT_SIZE 15
gmtset HEADER_OFFSET -1.5
gmtset ANNOT_FONT_SIZE 10

fMisL='Misfit2_L.out'
### plot Location Misfits ###
REG=-R0/1500/0.3/10
SCA=-JX10/4l
psbasemap -Ba500f100/a2f3:."Misfit**2 Location":WeSn $REG $SCA -X4 -Y15.5 -K > $psout
if [ -e $fMisL ]; then
	for iter in 1 2 3; do
		awk -v iter=$iter 'BEGIN{i=0}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $1,$2}}' $fMisL | psxy -R -J -A -W3,${color[$iter]} -O -K >> $psout
	done
fi

fMisF='Misfit2_F.out'
### plot Focal Misfits ###
REG=-R0/5000/0.3/300
SCA=-JX10/4l
psbasemap -Ba2000f500/a2f3:."Misfit**2 Focal":WeSn $REG $SCA -Y-5 -O -K >> $psout
if [ -e $fMisF ]; then
	for iter in 1 2 3; do
		awk -v iter=$iter 'BEGIN{i=0}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $1,$2}}' $fMisF | psxy -R -J -A -W1,${color[$iter]} -O -K >> $psout
	done
fi

perlst=( 10 16 22 30 40 )
iperbeg=0
ps=0.03
lw=2
### plot Group fit results ###
iter=$2
REG=-R0/360/-50/50
SCA=-JX10/8
psbasemap -Ba60f20/a20f5:."Fit Group":WeSn $REG $SCA -X12 -Y1 -O -K >> $psout
iper=$iperbeg
for per in ${perlst[@]}; do
   fsta=${type}_azi_data_pred_${per}sec.txt_sta
	fbin=${type}_azi_data_pred_${per}sec.txt_bin
   grep -v '\-12345' $fsta | awk 'NF>1' | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $4,$5-$6}}' | psxy -R -J -A -Sc${ps} -G${color[$iper]} -O -K >> $psout
   grep -v '\-12345' $fsta | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $4,$7}}' | psxy -R -J -A -W${lw},${color[$iper]} -O -K >> $psout
	grep -v '\-12345' $fbin | awk 'NF>1' | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $1,$2+$3,$4}}' | psxy -R -J -A -Ey0.2/${lw},${color[$iper]} -O -K >> $psout
   let iper++
done

### plot Phase fit results ###
REG=-R0/360/-30/30
SCA=-JX10/8
psbasemap -Ba60f20/a20f5:."Fit Phase":WeSn $REG $SCA -X-12 -Y-10 -O -K >> $psout
iper=$iperbeg
for per in ${perlst[@]}; do
   fsta=${type}_azi_data_pred_${per}sec.txt_sta
	fbin=${type}_azi_data_pred_${per}sec.txt_bin
   grep -v '\-12345' $fsta | awk 'NF>1' | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $4,$8-$9}}' | psxy -R -J -A -Sc${ps} -G${color[$iper]} -O -K >> $psout
   grep -v '\-12345' $fsta | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $4,$10}}' | psxy -R -J -A -W${lw},${color[$iper]} -O -K >> $psout
	grep -v '\-12345' $fbin | awk 'NF>1' | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $1,$5+$6,$7}}' | psxy -R -J -A -Ey0.2/${lw},${color[$iper]} -O -K >> $psout
   let iper++
done

### plot Ampplitude fit results ###
#REG=-R0/360/500/500000
REG=-R0/360/${ampmin}/${ampmax}
SCA=-JX10/8l
psbasemap -Ba60f20/a20f5:."Fit Amplitudes":WeSn $REG $SCA -X12 -O -K >> $psout
iper=$iperbeg
if [ $# == 5 ]; then
	ampfactor=$5
else
	ampfactor=1
fi
for per in ${perlst[@]}; do
   fsta=${type}_azi_data_pred_${per}sec.txt_sta
	fbin=${type}_azi_data_pred_${per}sec.txt_bin
   grep -v '\-12345' $fsta | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $4,$11}}' | psxy -R -J -A -Sc${ps} -G${color[$iper]} -O -K >> $psout
   grep -v '\-12345' $fsta | awk -v iter=$iter -v ampfactor=$ampfactor 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter&&NF>0){print $4,$12*ampfactor}}' | psxy -R -J -A -W${lw},${color[$iper]} -O -K >> $psout
	grep -v '\-12345' $fbin | awk 'NF>1' | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $1,$8+$9,$10}}' | psxy -R -J -A -Ey0.2/${lw},${color[$iper]} -O -K >> $psout
   let iper++
done


pwd | psxy -R -J -O >> $psout
echo $psout
