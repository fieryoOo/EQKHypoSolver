#!/bin/bash

CheckFile() {
	local fname=$1
	if [ ! -e $fname ]; then
		echo "Cannot find/access file "$fname
		exit
	fi
}

PlotText() {
	local _REG=$1
	local _text=$2
	local _XS=$3

	#local _lb=('a' 'b' 'c' 'd' 'e' 'f')
	#local _title=`echo ${_text} | awk -v lb=${_lb[ifile]} '{print "("lb")  "$0}'`
	local _title=$_text
	echo ${_title}
	local _llon=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v xs=$_XS '{print $1+(0.0+xs)*($2-$1)}'`
	local _ulat=`echo $_REG | sed s/'\-R'/''/ | awk -F/ '{print $4+0.02*($4-$3)}'`
	echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "15. 0. 20 LT", $0}' | pstext -R -J -Wlightgray,O3 -O -K -N >> $psout
	let ifile++
}

NormSac() {
	local _fsac=$1
	$exeSacnorm $_fsac $_fsac $tmin $tmax
}

PlotSac() {
	local _fsac=$1
	local _lflag=$2
	$exeSacdump $_fsac |	psxy -R -J $_lflag -O -K >> $psout
	#$exeSacnorm $_fsac $_fsac $tmin $tmax
	#rm -f ${_fsac}.txt
}


### main ###
if [ $# -lt 2 ]; then
	echo "Usage: "$0" [sac1] [sac2] ..."
	exit
fi

#sta=$1

# external programs
dircode=/home/tianye/code
if [ ! -e $dircode ]; then dircode=/projects/yeti4009/code; fi
if [ ! -e $dircode ]; then echo "code dir not found!"; exit; fi
dirscript=/home/tianye/EQKHypoSolver/Scripts
if [ ! -e $dirscript ]; then dirscript=/projects/yeti4009/eqkhyposolver/Scripts; fi
if [ ! -e $dirscript ]; then echo "script dir not found!"; exit; fi
exeFilter=`which SAC_filter`
exeSacdump=`which SAC_dump`
exeSacnorm=`which SAC_norm`
exeSacARatio=${dirscript}/SAC_ampratio
exePShift=${dircode}/Programs/WaveformFromPhaseDispersion/ComputePhaseShift
#exePssac=`which pssac`

# find sac files
fsacs=('' '')
for (( isac=1; isac<=$#; isac++ )); do
	fsac=${!isac}
	fsacs[isac]=$fsac
	CheckFile $fsac; SAC_touch $fsac
done
fsac1=${fsacs[1]}

#fsac1D=SAC_syn_1D/${sta}.BHZ_rev.SAC; CheckFile $fsac1D

# gmt setup
gmtset HEADER_FONT_SIZE 15
gmtset HEADER_OFFSET 0.
gmtset LABEL_FONT_SIZE 12
gmtset LABEL_OFFSET -0.15
gmtset ANNOT_FONT_SIZE 10
gmtset ANNOT_OFFSET 0.1

# line types
ltype_all=( '' -W8,100/100/100 -W5,steelblue,- -W5,lightbrown,- -W5,lightred,- -W3,red )
ltype=( '' ${ltype_all[1]} )
for (( isac=2; isac<=$#; isac++ )); do
	itype=`echo $# ${#ltype_all[@]} $isac | awk '{print $2-$1+$3-1}'`
	ltype[isac]=${ltype_all[itype]}
done

# plot sac files for comparisons
# broadband
# time range
dist=`saclst dist f $fsac1 | awk '{print $2}'`
tmin=`echo $dist | awk '{print $1/5.}'`
tmax=`echo $dist | awk '{print $1/2.}'`
if [ `echo $dist | awk '{if ($1<300){print 0}else{print 1}}'` == 0 ]; then
	tl=0
	tu=150
else
	tl=`echo $tmin | awk '{print $1}'`
	tu=`echo $tmax | awk '{print $1}'`
fi

# amp range
#Amaxold=`saclst depmin depmax f $fsac1 | awk '{max=-$2; if(max<$3)max=$3; print max}'`
#for (( isac=2; isac<=$#; isac++ )); do
#	fsac=${fsacs[isac]}
#	Amaxold=`saclst depmin depmax f $fsac | awk '{max=-$2; if(max<$3)max=$3; print max}' | awk -v Amaxold=$Amaxold '{if($1>Amaxold){print $1}else{print Amaxold}}'`
#done
#Amax=`echo $Amaxold | awk '{print $1*1.05}'`
#Atic=`echo $Amax | awk '{printf "%.0f", $1/2.5}'`
#REG=-R${tl}/${tu}/-${Amax}/${Amax}

# signal amplitude
mak=`saclst depmin depmax f ${fsac1} | awk '{max=$2;if(max<$3){max=$3}printf "%.2f", max}'`
amp=`echo $mak | awk '{printf "%.2f", $1*1.5}'`
tic=`echo $mak | awk '{printf "%.3f", $1*0.2}'`
REG=-R${tl}/${tu}/-${amp}/${amp}

SCA=-JX10/3.

# set up ps file
psout=${fsac1}.ps
pwd | psxy $REG $SCA -X6 -Y16.5 -K > $psout

psbasemap -R -J -Ba100f20/a${mak}f${tic}:."RealData(gray) - Synthectic(red)":Wesn -Y${YS} -O -K >> $psout
for (( isac=1; isac<=$#; isac++ )); do
	PlotSac ${fsacs[isac]} ${ltype[isac]}
done
PlotText $REG "broadband" 0.

# gauss filtered
# center periods to plot
perlst=( 8 12 18 )
YS=-3.; iper=1
for per in ${perlst[@]}; do
	# filter
	freq=`echo $per | awk '{print 1.0/$1}'`
	dfreq=`echo $per | awk '{e=2.718281828;per=$1; print 0.5*(1./e**( log(per)-0.15 ) - 1./e**( log(per)+0.15 ))}'`
	echo $freq $dfreq | awk '{print "current period(gaussian half range): "1./($1+$2), 1./$1, 1./($1-$2)}'
	for (( isac=1; isac<=$#; isac++ )); do
		fsac=${fsacs[isac]}
		$exeFilter $fsac -1. $freq $dfreq -1. ${fsac}_tmp
	done
	#$exeFilter $fsac1D -1. $freq $dfreq -1. ${fsac1D}_tmp
	# signal amplitude
	mak=`saclst depmin depmax f ${fsac1}_tmp | awk '{max=$2;if(max<$3){max=$3}printf "%.2f", max}'`
	amp=`echo $mak | awk '{printf "%.2f", $1*1.5}'`
	tic=`echo $mak | awk '{printf "%.3f", $1*0.2}'`
	REG=-R${tl}/${tu}/-${amp}/${amp}
	# plot sacs
	if [ $iper == 3 ]; then
		psbasemap $REG -J -Ba100f20:"Time (sec)":/a${mak}f${tic}:."":WeSn -Y${YS} -O -K >> $psout
	else
		psbasemap $REG -J -Ba100f20/a${mak}f${tic}Wesn -Y${YS} -O -K >> $psout
	fi
	for (( isac=1; isac<=$#; isac++ )); do
		#NormSac ${fsacs[isac]}_tmp
		PlotSac ${fsacs[isac]}_tmp ${ltype[isac]}
	done
	#$exePssac -R -J -W8,darkgray -O -K ${fsac1}_tmp -M1 >> $psout
	#$exePssac -R -J -W2,red -O -K ${fsac2}_tmp -M1 >> $psout
	# label
	PlotText $REG $per" sec" 0.
let iper++
done
echo ${fsacs[@]} | xargs -n1 | xargs -I file rm -f file'_tmp'

# plot amplitude transfer function
#REG=-R0.05/0.16/0/2.
#psbasemap $REG -JX8/5 -Ba0.05f0.01g0.05:"Frequency (Hz)":/a1f0.2g1:"Ratio"::."Spectral Ratio (Data/Syn)":WeSn -X11.8 -Y7. -O -K >> $psout
#PlotSac $fratio1D -W5,lightred,-
#PlotSac $fratio2 -W5,red
REG=-R6/20/0/2
psbasemap $REG -JX8/5 -Ba5f1g5:"Period (sec)":/a1f0.2g1:"Ratio"::."Spectral Ratio (Data/Syn)":WeSn -X11.8 -Y7. -O -K >> $psout
for (( isac=1; isac<=$#; isac++ )); do
	fsac=${fsacs[isac]}
	fratio=${fsac}_amp_ratio
	$exeSacARatio $fsac1 $fsac ${fratio} $tmin $tmax
	$exeSacdump $fratio | tac | awk '{if($1>0.)print 1./$1,$2}' | psxy -R -J ${ltype[isac]} -O -K >> $psout
	rm -f $fratio
done
sta=`echo $fsac1 | awk -F'/' '{print $NF}' | sed s/'.LH'/' '/ | sed s/'.BH'/' '/ | awk '{print $1}'`
PlotText $REG $sta 0.

# plot phase shift
REG=-R6/20/-5/5
psbasemap $REG -JX8/5 -Ba5f1g5:"Period (sec)":/a2f0.5g1:"Time (sec)"::."Phase Shift":WeSn -X0. -Y-7. -O -K >> $psout
# compute diff (file)
#fdiffs=('')
for (( isac=1; isac<=$#; isac++ )); do
	fsac=${fsacs[isac]}
	fdiff=${fsac}.diff
	$exePShift $fsac1 $fsac $fdiff
	psxy $fdiff -R -J -A ${ltype[isac]} -O -K >> $psout
done
PlotText $REG $sta 0.

# finalize
pwd | psxy -R -J -O >> $psout
echo $psout

