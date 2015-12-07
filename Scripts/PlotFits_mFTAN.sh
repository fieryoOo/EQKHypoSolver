#!/bin/bash

PlotText() {
   local _REG=$1
   local _text=$2
   local _XS=$3
   local _YS=$4

   #local _lb=('a' 'b' 'c' 'd' 'e' 'f')
   #local _title=`echo ${_text} | awk -v lb=${_lb[ifile]} '{print "("lb")  "$0}'`
   local _title=$_text
   echo ${_title}
   local _llon=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v xs=$_XS '{print $1+xs*($2-$1)}'`
   local _ulat=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v ys=$_YS '{print $4+ys*($4-$3)}'`
   echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "15. 0. 20 LT", $0}' | pstext -R -J -Wlightgray,O3 -O -K -N >> $psout
   #echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "12. 0. 20 LT", $0}' | pstext -R -J -O -K -N >> $psout
   let ifile++
}

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
psout=fits_${type}_mFTAN.ps 
color=(orange gold lightred steelblue midnightblue darkgray)
gmtset HEADER_FONT_SIZE 15
gmtset LABEL_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE 10
gmtset HEADER_OFFSET 0.
gmtset LABEL_OFFSET 0.
gmtset ANNOT_OFFSET 0.05

fMisL='Misfit2_L.out'
Emul=1
if [ $# -ge 5 ]; then Emul=$5; fi
### plot Location Misfits ###
REG=-R1/500/0.1/5
SCA=-JX10/4l
psbasemap -Ba100f20/a1f0.2WeSn $REG $SCA -X4 -Y15.5 -K > $psout
if [ -e $fMisL ]; then
	for iter in 1 2 3; do
		awk 'NF>0' $fMisL | awk -v iter=$iter 'BEGIN{flag=1;Nold=0}{if($2<Nold){flag++} if(flag==iter)print; Nold=$2}' | grep 'accepted' | awk -F"N = " '{print $2,$1}' | awk -v Emul=$Emul 'BEGIN{rcbest=99999}{is=$7;if(is<1){is=1} rc=$4*Emul/$1;print is,rc;if(rcbest>rc){rcbest=rc}}END{print 99999,rcbest}' | sort -g | psxy -R -J -A -W3,${color[$iter]} -O -K >> $psout
	done
fi
PlotText $REG "r-chi-square (Location)" 0.005 0

### legend for chi-square misfits
lw=6
#S dx1 symbol size fill pen [ dx2 text ]
pslegend -D500/5/3.4/2.1/TR -R -J -F -Ggray -O -K >> $psout <<- END
	G 0.1
	S 0.6 - 1. - ${lw},${color[1]} 1.3 iteration #1
	G 0.2
	S 0.6 - 1. - ${lw},${color[2]} 1.3 iteration #2
	G 0.2
	S 0.6 - 1. - ${lw},${color[3]} 1.3 iteration #3
END

fMisF='Misfit2_F.out'
### plot Focal Misfits ###
REG=-R1/10000/0.1/100
SCA=-JX10/4l
psbasemap -Ba2000f500/a5f1WeSn $REG $SCA -Y-5 -O -K >> $psout
if [ -e $fMisF ]; then
	for iter in 1 2 3; do
		#awk -v iter=$iter 'BEGIN{i=0}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $1,$2}}' $fMisF | psxy -R -J -A -W1,${color[$iter]} -O -K >> $psout
		awk 'NF>0' $fMisF | awk -v iter=$iter 'BEGIN{flag=1;Nold=0}{if($2<Nold){flag++} if(flag==iter)print; Nold=$2}' | grep 'accepted' | awk -F"N = " '{print $2,$1}' | awk -v Emul=$Emul 'BEGIN{rcbest=99999}{is=$7;if(is<1){is=1} rc=$4*Emul/$1;print is,rc;if(rcbest>rc){rcbest=rc}}END{print 99999,rcbest}' | sort -g | psxy -R -J -A -W3,${color[$iter]} -O -K >> $psout
	done
fi
PlotText $REG "r-chi-square (Focal)" 0.005 0

#pwd | psxy -R -J -O >> $psout


### fits of group, phase, and amplitude
perlst=( 8 10 16 22 30 40 )
fsource=${type}_source_patterns.txt
iperbeg=0
ps=0.05
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
	Niterlast=`grep '#' $fsta | wc -l | awk '{print $1-1}'`
   grep -v '\-12345' $fsta | awk 'NF>1' | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $4,$5-$6}}' | psxy -R -J -A -Sc${ps} -G${color[$iper]} -O -K >> $psout
	if [ $iter == $Niterlast ] && [ -e $fsource ]; then
		awk -v per=$per '$5==per{print $1,$2}' $fsource | psxy -R -J -A -W${lw},${color[$iper]} -O -K >> $psout
	else
		grep -v '\-12345' $fsta | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $4,$7}}' | psxy -R -J -A -W${lw},${color[$iper]} -O -K >> $psout
	fi
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
	Niterlast=`grep '#' $fsta | wc -l | awk '{print $1-1}'`
   grep -v '\-12345' $fsta | awk 'NF>1' | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $4,$8-$9}}' | psxy -R -J -A -Sc${ps} -G${color[$iper]} -O -K >> $psout
	if [ $iter == $Niterlast ] && [ -e $fsource ]; then
		awk -v per=$per '$5==per{print $1,$3}' $fsource | psxy -R -J -A -W${lw},${color[$iper]} -O -K >> $psout
	else
		grep -v '\-12345' $fsta | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $4,$10}}' | psxy -R -J -A -W${lw},${color[$iper]} -O -K >> $psout
	fi
	grep -v '\-12345' $fbin | awk 'NF>1' | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $1,$5+$6,$7}}' | psxy -R -J -A -Ey0.2/${lw},${color[$iper]} -O -K >> $psout
   let iper++
done

### plot Ampplitude fit results ###
#REG=-R0/360/500/500000
REG=-R0/360/${ampmin}/${ampmax}
SCA=-JX10/8l
psbasemap -Ba60f20/a2f3:."Fit Amplitudes":WeSn $REG $SCA -X12 -O -K >> $psout
#SCA=-JX10/8
#psbasemap -Ba60f20/a10000f5000:."Fit Amplitudes":WeSn $REG $SCA -X12 -O -K >> $psout
iper=$iperbeg
if [ $# -ge 6 ]; then
	ampfactor=$6
else
	ampfactor=1
fi
for per in ${perlst[@]}; do
   fsta=${type}_azi_data_pred_${per}sec.txt_sta
	fbin=${type}_azi_data_pred_${per}sec.txt_bin
	Niterlast=`grep '#' $fsta | wc -l | awk '{print $1-1}'`
   grep -v '\-12345' $fsta | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $4,$11}}' | psxy -R -J -A -Sc${ps} -G${color[$iper]} -O -K >> $psout
	if [ 0 == 1 ] && [ $iter == $Niterlast ] && [ -e $fsource ]; then
		awk -v per=$per -v ampfactor=$ampfactor '$5==per{print $1,$4*ampfactor}' $fsource | psxy -R -J -A -W${lw},${color[$iper]} -O -K >> $psout
	else
		grep -v '\-12345' $fsta | awk -v iter=$iter -v ampfactor=$ampfactor 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter&&NF>0){print $4,$12*ampfactor}}' | psxy -R -J -A -W${lw},${color[$iper]} -O -K >> $psout
	fi
	grep -v '\-12345' $fbin | awk 'NF>1' | awk -v iter=$iter 'BEGIN{i=-1}{if(substr($1,0,1)=="#"){i++}else if(i==iter){print $1,$8*(1.0+$9),$8*$10}}' | psxy -R -J -A -Ey0.2/${lw},${color[$iper]} -O -K >> $psout
   let iper++
done

pwd | psxy -R -J -O >> $psout
echo $psout

