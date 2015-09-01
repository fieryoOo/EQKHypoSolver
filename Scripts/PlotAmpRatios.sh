#!/bin/bash

GetRatio() {
	local _sac1=$1
	local _sac2=$2
	local _tmin=$3
	local _tmax=$4
	local _per=$5
	exeSAR=/home/tianye/EQKHypoSolver/Scripts/SAC_ampratio

	#Usage: ../../Scripts/SAC_ampratio [sac1] [sac2] [sac_outname (cout=print to screen)] [tmin(optional)] [tmax(optional)]
	$exeSAR $_sac1 $_sac2 cout $_tmin $_tmax | awk -v per=$_per '{f=1./per; if($1>f){print rl+($2-rl)*(f-fl)/($1-fl); exit} fl=$1;rl=$2}'
}


### main ###
if [ $# != 1 ]; then
	echo "Usage: "$0" [ SAC dir ]"
	exit
fi

### parameters
dir=$1; dir=${dir%/}
per=10
fcpt=/home/tianye/EQKHypoSolver/Scripts/AmpRatio.cpt

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

# plot, respectively, for the Z and T components
for comp in Z T; do
	pwd | psbasemap -R -J -Ba5f1:."Amp ratio ($comp component)":WeSn -O -K -X12 >> $psout
	for sacf in `ls ${dir}/*.LH${comp}_real.SAC`; do
		sacfs=`echo $sacf | sed s/'real'/'syn'/`
		if [ ! -e $sacfs ]; then
			echo "Warning: synthetic sac "$sacfs" not found!!! skipped."
			continue;
		fi
		loc=`$exeSPH $sacf stlo stla`
		trange=`$exeSPH $sacf dist e | awk '{tmin=$1/4.5-50.; if(tmin<0.){tmin=0.} tmax=$1/2.5+50.; if(tmax>$2){tmax=$2} print tmin, tmax}'`
		ratio=`GetRatio $sacf $sacfs $trange $per`
		echo $loc $ratio | psxy -R -J -Sc0.5 -W1 -C$fcpt -O -K >> $psout
	done # sacf
	psscale  -C$fcpt -B:"${_label}": -L -D5./-1./10./0.4h -O -K >> $psout
done # comp

pwd | psxy -R -J -O >> $psout
echo $psout
