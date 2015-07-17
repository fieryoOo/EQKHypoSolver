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
if [ $# != 2 ]; then
	echo "Usage: "$0" [sac1] [sac2]"
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
fsac1=$1; CheckFile $fsac1
fsac2=$2; CheckFile $fsac2
#fsac1D=SAC_syn_1D/${sta}.BHZ_rev.SAC; CheckFile $fsac1D

# compute diff (file)
fdiff2=${fsac2}.diff
$exePShift $fsac1 $fsac2 $fdiff2
#fdiff1D=${fsac1D}.diff

# gmt setup
gmtset HEADER_FONT_SIZE 15
gmtset HEADER_OFFSET 0.
gmtset LABEL_FONT_SIZE 12
gmtset LABEL_OFFSET -0.15
gmtset ANNOT_FONT_SIZE 10
gmtset ANNOT_OFFSET 0.1

# center periods to plot
perlst=( 8 12 18 )
# set up ps file
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
Amax1=`saclst depmin depmax f $fsac1 | awk '{max=-$2; if(max<$3)max=$3; print max}'`
Amax=`saclst depmin depmax f $fsac2 | awk '{max=-$2; if(max<$3)max=$3; print max}' | awk -v Amax1=$Amax1 '{if($1>Amax1){print $1}else{print Amax1}}' | awk '{print $1*1.05}'`

REG=-R${tl}/${tu}/-${Amax}/${Amax}
SCA=-JX10/3.

psout=${fsac1}.ps
pwd | psxy $REG $SCA -X6 -Y16.5 -K > $psout

# plot sac files for comparisons
YS=-3.
iper=1

# normalize both sacs
#NormSac ${fsac1}
#NormSac ${fsac2}
#NormSac ${fsac1D}

# broadband
ltype=( -W8,100/100/100 -W5,lightred,- -W3,red )
Atic=`echo $Amax | awk '{printf "%.0f", $1/2.5}'`
psbasemap -R -J -Ba50f10/a${Atic}f${Atic}:."RealData(gray) - Synthectic(red)":Wesn -Y${YS} -O -K >> $psout
PlotSac ${fsac1} ${ltype[0]}
#PlotSac ${fsac1D} ${ltype[1]}
PlotSac ${fsac2} ${ltype[2]}
PlotText $REG "broadband" 0.
# gauss filtered
for per in ${perlst[@]}; do
	# filter
	freq=`echo $per | awk '{print 1.0/$1}'`
	dfreq=`echo $per | awk '{e=2.718281828;per=$1; print 0.5*(1./e**( log(per)-0.15 ) - 1./e**( log(per)+0.15 ))}'`
	echo $freq $dfreq | awk '{print "current period(gaussian half range): "1./($1+$2), 1./$1, 1./($1-$2)}'
	$exeFilter $fsac1 -1. $freq $dfreq -1. ${fsac1}_tmp
	$exeFilter $fsac2 -1. $freq $dfreq -1. ${fsac2}_tmp
	#$exeFilter $fsac1D -1. $freq $dfreq -1. ${fsac1D}_tmp
	# signal amplitude
	mak=`saclst depmin depmax f ${fsac2}_tmp | awk '{max=$2;if(max<$3){max=$3}printf "%.2f", max}'`
	amp=`echo $mak | awk '{printf "%.2f", $1*1.5}'`
	tic=`echo $mak | awk '{printf "%.3f", $1*0.2}'`
	REG=-R${tl}/${tu}/-${amp}/${amp}
	# plot sacs
	if [ $iper == 3 ]; then
		psbasemap $REG -J -Ba100f20:"Time (sec)":/a${mak}f${tic}:."":WeSn -Y${YS} -O -K >> $psout
	else
		psbasemap $REG -J -Ba100f20/a${mak}f${tic}Wesn -Y${YS} -O -K >> $psout
	fi
	PlotSac ${fsac1}_tmp ${ltype[0]}
	#PlotSac ${fsac1D}_tmp ${ltype[1]}
	PlotSac ${fsac2}_tmp ${ltype[2]}
	#$exePssac -R -J -W8,darkgray -O -K ${fsac1}_tmp -M1 >> $psout
	#$exePssac -R -J -W2,red -O -K ${fsac2}_tmp -M1 >> $psout
	# label
	PlotText $REG $per" sec" 0.
let iper++
done
rm -f ${fsac1}_tmp ${fsac2}_tmp ${fsac1D}_tmp

# plot amplitude transfer function
fratio2=${fsac2}_amp_ratio
#fratio1D=${fsac1D}_amp_ratio
$exeSacARatio $fsac1 $fsac2 ${fratio2} $tmin $tmax
#$exeSacARatio $fsac1 $fsac1D ${fratio1D} $tmin $tmax
#REG=-R0.05/0.16/0/2.
#psbasemap $REG -JX8/5 -Ba0.05f0.01g0.05:"Frequency (Hz)":/a1f0.2g1:"Ratio"::."Spectral Ratio (Data/Syn)":WeSn -X11.8 -Y7. -O -K >> $psout
#PlotSac $fratio1D -W5,lightred,-
#PlotSac $fratio2 -W5,red
REG=-R6/20/0/2
psbasemap $REG -JX8/5 -Ba5f1g5:"Period (sec)":/a1f0.2g1:"Ratio"::."Spectral Ratio (Data/Syn)":WeSn -X11.8 -Y7. -O -K >> $psout
$exeSacdump $fratio2 | tac | awk '{if($1>0.)print 1./$1,$2}' | psxy -R -J -W5,red -O -K >> $psout
sta=`echo $fsac1 | awk -F'/' '{print $NF}' | sed s/'.LH'/' '/ | sed s/'.BH'/' '/ | awk '{print $1}'`
PlotText $REG $sta 0.
rm -f $fratio2 $fratio1D

# plot phase shift
REG=-R6/20/-5/5
psbasemap $REG -JX8/5 -Ba5f1g5:"Period (sec)":/a2f0.5g1:"Time (sec)"::."Phase Shift":WeSn -X0. -Y-7. -O -K >> $psout
#psxy $fdiff1D -R -J -A -W8,lightblue,- -O -K >> $psout
psxy $fdiff2 -R -J -A -W8,steelblue -O -K >> $psout
PlotText $REG $sta 0.

# finalize
pwd | psxy -R -J -O >> $psout
echo $psout

