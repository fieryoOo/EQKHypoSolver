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
   echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "10. 0. 20 LT", $0}' | pstext -R -J -Wlightgray,O3 -O -K -N >> $psout
   #echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "12. 0. 20 LT", $0}' | pstext -R -J -O -K -N >> $psout
   #let ifile++
}

PlotSingle() {
	local _fin=$1
	local _type=$2 # 1=group, 2=phase, 3=amplitude
	local _nos=$3 # 1=solution1, 2=solution2, 3=solution3
	local _XS=$4
	local _YS=$5
	local _col; local _REG; local _mrk;
	local _lt=("" "-W3,blue" "-W3,darkgray" "-W3,red")
	if [ $_type == 1 ]; then
		_col=2; _REG=-R0/360/-${per}/${per}; _mrk=$per; _typename='Group'
	elif [ $_type == 2 ]; then
		_col=3; _REG=-R0/360/-${per}/${per}; _mrk=$per; _typename='Phase'
	elif [ $_type == 3 ]; then
		_col=4; _REG=-R0/360/0/0.2; _mrk=0.2; _typename='Amplitude'
	else
		echo "Wrong type! "$_type
		exit
	fi
	local _SCA=-JX5/3.7
	local _tic=`echo $_mrk | awk '{print $1*0.2}'`
	if [ $_nos == 1 ]; then 
		psbasemap $_REG $_SCA -Ba120f20/a${_mrk}f${_tic}WeSn -X$_XS -Y$_YS -O -K >> $psout; 
		PlotText $_REG "$per sec $_typename" 0 0.05
	fi
	awk -v col=$_col '{print $1,$(col)}' $_fin | psxy -R -J ${_lt[_nos]} -O -K >> $psout
}

### main ###
gmtset ANNOT_FONT_SIZE 8
gmtset ANNOT_OFFSET 0.05

fsp1=solution_blue; ftmp1=${fsp1}.tmp
fsp2=solution_gray; ftmp2=${fsp2}.tmp
fsp3=solution_red ; ftmp3=${fsp3}.tmp

psout=SourcePatternComparison.ps
pwd | psxy -Rg -JX1 -X14 -Y29 -K -P > $psout

for per in 5 8 10 16 25 40; do
	awk -v per=$per '$5==per' $fsp1 > $ftmp1
	awk -v per=$per '$5==per' $fsp2 > $ftmp2
	awk -v per=$per '$5==per' $fsp3 > $ftmp3

	PlotSingle $ftmp1 1 1 -12 -4.4
	PlotSingle $ftmp2 1 2
	PlotSingle $ftmp3 1 3

	PlotSingle $ftmp1 2 1 6 0
	PlotSingle $ftmp2 2 2
	PlotSingle $ftmp3 2 3

	PlotSingle $ftmp1 3 1 6 0
	PlotSingle $ftmp2 3 2
	PlotSingle $ftmp3 3 3
done

pwd | psxy -R -J -O >> $psout
rm -f $ftmp1 $ftmp2 $ftmp3
echo $psout
