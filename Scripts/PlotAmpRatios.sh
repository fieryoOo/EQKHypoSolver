#!/bin/bash

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

GetRatio() {
	local _sac1=$1
	local _sac2=$2
	local _tmin=$3
	local _tmax=$4
	local _per=$5
	exeSAR=${dirS}/SAC_ampratio

	#Usage: ../../Scripts/SAC_ampratio [sac1] [sac2] [sac_outname (cout=print to screen)] [tmin(optional)] [tmax(optional)]
	$exeSAR $_sac1 $_sac2 cout $_tmin $_tmax | awk -v per=$_per '{f=1./per; if($1>f){print rl+($2-rl)*(f-fl)/($1-fl); exit} fl=$1;rl=$2}'
}


### main ###
if [ $# != 2 ] && [ $# != 3 ]; then
	echo "Usage: "$0" [waveform dir] [per] [fcpt (optional)]"
	exit
fi

### parameters
dir=$1; dir=${dir%/}
per=$2
dirS=/home/tianye/EQKHypoSolver/Scripts
if [ ! -e $dirS ]; then
	dirS=/projects/yeti4009/eqkhyposolver/Scripts
fi
fcpt=${dirS}/AmpRatio.cpt
if [ $# == 3 ]; then fcpt=$3; fi

### gmt settings
gmtset HEADER_FONT_SIZE 15
gmtset LABEL_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE 10
gmtset HEADER_OFFSET 0
gmtset LABEL_OFFSET 0
gmtset ANNOT_OFFSET 0.05

### output ps file
psout=${dir}_AmpRatios.ps; rm -f $psout
REG=-R-123/-108/33/46
SCA=-Jn-115/0.7
pwd | psxy $REG $SCA -K -X-7 -Y5 > $psout

### exe s
exeSPH=/home/tianye/code/Programs/SACOperation/SAC_printHD
if [ ! -e $exeSPH ]; then
	exeSPH=/projects/yeti4009/code/Programs/SACOperation/SAC_printHD
fi

# plot, respectively, for the Z and T components
for comp in Z T; do
	sum_ratio=0; Nsac=0
	pwd | psbasemap -R -J -Ba5f1:."Amp ratio (obs/syn $comp comp)":WeSn -O -K -X12 >> $psout
	for sacf in `ls ${dir}/*.LH${comp}_real.SAC`; do
		sacfs=`echo $sacf | sed s/'real'/'syn'/`
		if [ ! -e $sacfs ]; then
			echo "Warning: synthetic sac "$sacfs" not found!!! skipped."
			continue;
		fi
		loc=`$exeSPH $sacf stlo stla`
		trange=`$exeSPH $sacf dist e | awk '{tmin=$1/4.5-50.; if(tmin<0.){tmin=0.} tmax=$1/2.5+50.; if(tmax>$2){tmax=$2} print tmin, tmax}'`
		ratio=`GetRatio $sacf $sacfs $trange $per`
		sum_ratio=`echo $sum_ratio $ratio | awk '{print $1+$2}'`; let Nsac++;
echo $loc $ratio
		echo $loc $ratio | psxy -R -J -Sc0.5 -W1 -C$fcpt -O -K >> $psout
	done # sacf
	psscale  -C$fcpt -B:"${_label}": -L -D5./-1./10./0.4h -O -K >> $psout
	avg_ratio=`echo $sum_ratio $Nsac | awk '{print $1/$2}'`
	PlotText $REG "average: $avg_ratio" 0.01
done # comp

pwd | psxy -R -J -O >> $psout
echo $psout
