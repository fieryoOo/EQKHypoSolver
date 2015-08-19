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
	local _YS=$4

	#local _lb=('a' 'b' 'c' 'd' 'e' 'f')
	#local _title=`echo ${_text} | awk -v lb=${_lb[ifile]} '{print "("lb")  "$0}'`
	local _title=$_text
	echo ${_title}
	#local _llon=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v xs=$_XS '{print $1+0.0*($2-$1)+xs}'`
	#local _ulat=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v ys=$_YS '{print $4+0.0*($4-$3)+ys}'`
	local _llon=$_XS
	local _ulat=$_YS
	#echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "15. 0. 20 LT", $0}' | pstext -R -J -Wwhite,O3 -O -K -N >> $psout
	echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "10. 0. 20 LT", $0}' | pstext -R -J -Wwhite,O3 -O -K -N >> $psout
	let ifile++
}

PlotSac() {
	local _fsac=$1
	local _dis=$2
	local _lflag=$3
	$exeSacdump $_fsac | awk -v dis=$_dis '{print $1,$2*0.005+dis}' |	psxy -R -J $_lflag -O -K >> $psout
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
#REG=-R0/400/-15000/15000
REG=-R0/450/-100/1100
#REG=-R0/450/300/700
SCA=-JX15/20.
pwd | psxy $REG $SCA -X3.5 -Y5 -K -P > $psout

# sort flst
ftmp=${flst}.PlotSacs.tmp
rm -f $ftmp
while read fsacr fsacs1 fsacs2 fsacs3 stmp; do
	sta=`echo $fsacs1 | awk -F'/' '{print $NF}' | sed s/'.LH'/' '/ | sed s/'.BH'/' '/ | awk '{print $1}'`
	SAC_touch $fsacs1
	dis=`saclst dist f $fsacs1 | awk '{print $2}'`
	echo $sta $dis $fsacr $fsacs1 $fsacs2 $fsacs3 >> $ftmp
done < $flst
sort -g -k2 $ftmp > ${ftmp}.tmp
mv ${ftmp}.tmp $ftmp

# gmt setting
gmtset HEADER_FONT_SIZE 18
gmtset LABEL_FONT_SIZE 15
gmtset ANNOT_FONT_SIZE 12
gmtset HEADER_OFFSET 0
gmtset LABEL_OFFSET 0
gmtset ANNOT_OFFSET 0.1

# plot sac files
#-W8,100/100/100 -W5,lightred,- -W3,red
pwd | psbasemap -R -J -Ba200f40:"time (sec)":/a500f100:"distance (km)"::."Waveforms (Real-3D-SLU)":WeSn -O -K >> $psout
colors=( '' 100/100/100 red 60/200/30 lightbrown lightred )
ltypes=( '' -W5,${colors[1]} -W3,${colors[2]} -W2,${colors[3]} -W2,${colors[4]},- -W3,${colors[5]},- )
while read sta dis fsacr fsacs1 fsacs2 fsacs3; do
	#pwd | psbasemap -R -J -Ba200f40/a15000f3000WeSn -Y-2. -O -K >> $psout
	CheckFile $fsacr; PlotSac $fsacr $dis ${ltypes[1]}
	CheckFile $fsacs1; PlotSac $fsacs1 $dis ${ltypes[2]}
	CheckFile $fsacs3; PlotSac $fsacs3 $dis ${ltypes[4]}
	CheckFile $fsacs2; PlotSac $fsacs2 $dis ${ltypes[3]}
	#tu=`echo $dis | awk '{print int($1*0.5+0.5)+25}'`
	#tl=`echo $dis | awk '{print int($1*0.2+0.5)-25}'`
	tl=`SAC_printHD $fsacs1 b | awk '{print $1-30}'`
	tu=`SAC_printHD $fsacs1 e`
	if [ $tu -gt 150 ]; then
		PlotText $REG $sta ${tl} $dis
	else
		PlotText $REG $sta ${tu} $dis
	fi
done < $ftmp

# wave color legend
#S dx1 symbol size fill pen [ dx2 text ]
#S 0.2i - 0.3i - 2p,${colors[1]},4:0p 0.5i Real
pslegend -R -D-3/970/2.8/2.4/BL -J -F -G220 -O -K <<- EOF >> $psout
G 0.05i
S 0.25i - 0.4i - 1.5p,${colors[1]} 0.6i Real
G 0.1i
S 0.25i - 0.4i - 1p,${colors[2]} 0.6i 3D
G 0.1i
S 0.25i - 0.4i - 1p,${colors[3]} 0.6i SLU
G 0.05i
EOF

# plot ends
rm -f $ftmp
pwd | psxy -R -J -O >> $psout
echo $psout

