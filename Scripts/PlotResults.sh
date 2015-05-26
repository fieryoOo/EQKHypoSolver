#!/bin/bash

### main ###
if [ $# != 4 ] && [ $# != 5 ] && [ $# != 6 ]; then
   echo "Usage: "$0" [type(R/L)] [iter#] [minamp] [maxamp] [E_multiplier (optional)] [amp factor (optional)]"
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
color=(orange lightred forestgreen steelblue midnightblue darkgray)
gmtset HEADER_FONT_SIZE 15
gmtset HEADER_OFFSET -1.5
gmtset ANNOT_FONT_SIZE 10

fMisL='Misfit2_L.out'
Emul=1
if [ $# -ge 5 ]; then Emul=$5; fi
### plot Location Misfits ###
REG=-R1/500/0.1/5
SCA=-JX10l/4
psbasemap -Ba100f20/a1f0.2:."Misfit**2 Location":WeSn $REG $SCA -X4 -Y15.5 -K > $psout
if [ -e $fMisL ]; then
	for iter in 1 2 3; do
		awk 'NF>0' $fMisL | awk -v iter=$iter 'BEGIN{flag=1;Nold=0}{if($2<Nold){flag++} if(flag==iter)print; Nold=$2}' | grep 'accepted' | awk -F"N = " '{print $2,$1}' | awk -v Emul=$Emul 'BEGIN{rcbest=99999}{is=$7;if(is<1){is=1} rc=$4*Emul/$1;print is,rc;if(rcbest>rc){rcbest=rc}}END{print 99999,rcbest}' | sort -g | psxy -R -J -A -W3,${color[$iter]} -O -K >> $psout
	done
fi

fMisF='Misfit2_F.out'
### plot Focal Misfits ###
REG=-R1/10000/0.1/20
SCA=-JX10/4
psbasemap -Ba2000f500/a5f1:."Misfit**2 Focal":WeSn $REG $SCA -Y-5 -O -K >> $psout
if [ -e $fMisF ]; then
	for iter in 1 2 3; do
		#awk -v iter=$iter 'BEGIN{i=0}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $1,$2}}' $fMisF | psxy -R -J -A -W1,${color[$iter]} -O -K >> $psout
		awk 'NF>0' $fMisF | awk -v iter=$iter 'BEGIN{flag=1;Nold=0}{if($2<Nold){flag++} if(flag==iter)print; Nold=$2}' | grep 'accepted' | awk -F"N = " '{print $2,$1}' | awk -v Emul=$Emul 'BEGIN{rcbest=99999}{is=$7;if(is<1){is=1} rc=$4*Emul/$1;print is,rc;if(rcbest>rc){rcbest=rc}}END{print 99999,rcbest}' | sort -g | psxy -R -J -A -W3,${color[$iter]} -O -K >> $psout
	done
fi

perlst=( 8 10 16 22 30 40 )
iperbeg=0
ps=0.1
lw=2
### plot Group fit results ###
iter=$2
#REG=-R0/360/-50/50
REG=-R0/360/-25/25
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
#REG=-R0/360/-30/30
REG=-R0/360/-15/15
SCA=-JX10/8
psbasemap -Ba60f20/a15f5:."Fit Phase":WeSn $REG $SCA -X-12 -Y-10 -O -K >> $psout
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
if [ $# -ge 6 ]; then
	ampfactor=$6
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
