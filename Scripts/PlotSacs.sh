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

PlotSac() {
	local _fsac=$1
	local _lflag=$2
	$exeSacdump $_fsac |	psxy -R -J $_lflag -O -K >> $psout
	#$exeSacnorm $_fsac $_fsac $tmin $tmax
	#rm -f ${_fsac}.txt
}


### main ###
if [ $# != 1 ]; then
	echo "Usage: "$0" [saclist]"
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

# plot starts
flst=$1
psout=${flst}.ps
REG=-R0/400/-15000/15000
SCA=-JX10/2.
pwd | psxy $REG $SCA -X3.5 -Y28 -K -P > $psout

# sort flst
ftmp=${flst}.PlotSacs.tmp
rm -f $ftmp
while read fsacr fsacs stmp; do
	sta=`echo $fsacs | awk -F'/' '{print $NF}' | sed s/'.LH'/' '/ | sed s/'.BH'/' '/ | awk '{print $1}'`
	SAC_touch $fsacs
	dis=`saclst dist f $fsacs | awk '{print $2}'`
	echo $dis $fsacr $fsacs $sta >> $ftmp
done < $flst
sort -g -k1 $ftmp | awk '{print $2,$3,$4}' > ${ftmp}.tmp
mv ${ftmp}.tmp $ftmp

# plot sac files
#-W8,100/100/100 -W5,lightred,- -W3,red
while read fsacr fsacs sta stmp; do
	pwd | psbasemap -R -J -Ba200f40/a15000f3000WeSn -Y-2. -O -K >> $psout
	CheckFile $fsacr
	PlotSac $fsacr -W5,100/100/100
	CheckFile $fsacs
	PlotSac $fsacs -W2,red
	PlotText $REG $sta 0.
done < $ftmp

# plot ends
rm -f $ftmp
pwd | psxy -R -J -O >> $psout
echo $psout

