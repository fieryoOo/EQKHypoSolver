#!/bin/bash


PlotText() {
	local _REG=$1
	local _text=$2
	local _XS=$3

	#local _lb=('a' 'b' 'c' 'd' 'e' 'f')
	#local _title=`echo ${_text} | awk -v lb=${_lb[ifile]} '{print "("lb")  "$0}'`
	local _title=$_text
	echo ${_title}
	local _llon=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v xs=$_XS '{print $1+(0.01+xs)*($2-$1)}'`
	local _ulat=`echo $_REG | sed s/'\-R'/''/ | awk -F/ '{print $4+0.02*($4-$3)}'`
	echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "15. 0. 20 LT", $0}' | pstext -R -J -Wlightgray,O3 -O -K -N >> $psout
	let ifile++
}

ComputeMeanStd() {
   h_infile=$1
   columnN=$2
   mean=`awk -v cN=$columnN 'BEGIN{sum=0}{sum+=$cN}END{print sum/NR}' $h_infile`
   std=`awk -v cN=$columnN -v mean=$mean 'BEGIN{sum=0}{sum+=($cN-mean)**2}END{print sqrt(sum/(NR-1))}' $h_infile`
	echo -n $mean" "$std
}


GetRegion() {
	local ArrayName=$1
	local _lonmin=360.
	local _lonmax=-360.
	local _latmin=90.
	local _latmax=-90.
	# find min max
	for loc in `echo ${!ArrayName}`; do
		lon=`echo $loc | awk -F_ '{print $1}'`
		lat=`echo $loc | awk -F_ '{print $2}'`
		_lonmin=`echo $lon $_lonmin | awk '{print $1<$2 ? $1 : $2}'`
		_lonmax=`echo $lon $_lonmax | awk '{print $2<$1 ? $1 : $2}'`
		_latmin=`echo $lat $_latmin | awk '{print $1<$2 ? $1 : $2}'`
		_latmax=`echo $lat $_latmax | awk '{print $2<$1 ? $1 : $2}'`
	done
	# extend by 20% on each side
	local _lon1=`echo $_lonmin $_lonmax | awk '{print $1-($2-$1)*0.4}'`
	local _lon2=`echo $_lonmin $_lonmax | awk '{print $2+($2-$1)*0.4}'`
	local _lat1=`echo $_latmin $_latmax | awk '{print $1-($2-$1)*0.4}'`
	local _lat2=`echo $_latmin $_latmax | awk '{print $2+($2-$1)*0.4}'`
	echo -R${_lon1}/${_lon2}/${_lat1}/${_lat2}
}


ExtendRegion() {
	# input
	local _lon=`echo $1 | awk '{lon=$1; if(lon<0.){lon+=360.} print lon}'`
	local _lat=$2
	local _stdlon=0
	local _stdlat=0
	if [ $# == 4 ]; then
		_stdlon=$3
		_stdlat=$4
	fi
	# current
	local _lonl=`echo $REG | cut -d/ -f1 | sed s/'-R'/''/`
	local _lonu=`echo $REG | cut -d/ -f2`
	local _latl=`echo $REG | cut -d/ -f3`
	local _latu=`echo $REG | cut -d/ -f4`
	# update (extended by ~2 km on each direction)
	_lonl=`echo $_lonl $_lon $_stdlon | awk '{lon2=$2-2.*$3-0.01; if($1<lon2){print $1}else{print lon2}}'`
	_lonu=`echo $_lonu $_lon $_stdlon | awk '{lon2=$2+2.*$3+0.01; if($1>lon2){print $1}else{print lon2}}'`
	_latl=`echo $_latl $_lat $_stdlat | awk '{lat2=$2-2.*$3-0.01; if($1<lat2){print $1}else{print lat2}}'`
	_latu=`echo $_latu $_lat $_stdlat | awk '{lat2=$2+2.*$3+0.01; if($1>lat2){print $1}else{print lat2}}'`
	REG=-R${_lonl}/${_lonu}/${_latl}/${_latu}
}


PlotEllipse() {
	local _psout=$1
	local _clon=`echo $2 | awk '{lon=$1; if(lon<0.){lon+=360.} print lon}'`
	local _clat=$3
	local DistEXE=`which get_dist`
	local _lonfactor=`$DistEXE $_clat 0 $_clat 1 d`
	local _latfactor=`$DistEXE 0 $_clon 1 $_clon d`
	# diameters of the ellipse ( +/- 2*std-dev )
	local _Dlon=`echo $4 $_lonfactor | awk '{print $1*$2*4.}'`
	local _Dlat=`echo $5 $_latfactor | awk '{print $1*$2*4.}'`
	local _rclon=$6
	local _rclat=$7
	local _color=$8
	echo -e $_rclon $_rclat"\n"$_clon $_clat | psxy -R -J -W6,${_color} -O -K >> $psout
	echo $_clon $_clat 90. $_Dlon $_Dlat | psxy -R -J -SE -W6,${_color} -O -K >> $psout
}



PlotPoint() {
	local _psout=$1
	local _lon=`echo $2 | awk '{lon=$1; if(lon<0.){lon+=360.} print lon}'`
	local _lat=$3
	local _colorFill=$4
	local _colorPen=$5
	echo $_lon $_lat | psxy -R -J -Sa.5 -G$_colorFill -W2,$_colorPen -O -K >> $_psout
}


PlotBeachball() { 
   local _psout=$1 
   local _stk=$2
   local _dip=$3
   local _rak=$4
   local _dep=`echo $5 | awk '{printf "%.1f",$1}'`  # rounded depth
	local _lon=`echo $6 | awk '{lon=$1; if(lon<0.){lon+=360.} print lon}'`
   local _lat=$7
   local _colorFill=$8
   local _colorPen=$9
	local _text=${10}

	local _bsize=0.6
	if [ `echo $_stk $_dip $_rak | awk '{if($1==0&&$2==0&&$3==0){print 0}else{print 1}}'` == 1 ]; then
		echo $_lon $_lat $_dep $_stk $_dip $_rak 5. 0 0 | psmeca -R -J -Sa${_bsize} -G$_colorFill -W2,$_colorPen -O -K >> $_psout
	else
		echo $_lon $_lat | psxy -R -J -Sc${_bsize} -G$_colorFill -W2,$_colorPen -O -K >> $_psout
	fi
	if [ $# -gt 9 ]; then
		echo $_lon $_lat $_text | awk '{print $1, $2-0.002, "15. 0. 20 LT", $3}' | pstext -R -J -Wlightgray,O3 -O -K -N >> $psout
	fi
}


### main ###
#if [ $# != 1 ]; then
#	echo "Usage: "$0" [event name]"
#	exit
#fi

#ename=`echo $1 | cut -d/ -f1`
ename=`pwd | awk -F/ '{print $NF}'`

# region
REG=-R360/-360/90/-90

idir=0
lonavg=0
for rdir in `ls -d results_SAMC_?_1000_Ei results_SAMC_?_1000_1D`; do
#for rdir in `ls -d results_SAMC_?_1000_Ei`; do

	type[idir]=`echo $rdir | awk -F/ '{print $NF}' | cut -d_ -f3,4,5`

	### location range from MC (PosteriorD.txt)
	fMC=${rdir}/PosteriorD.txt
	if [ ! -e $fMC ]; then
		echo "   Warning: file "$fMC" not found"
		continue
	fi
	# extract data: discard the initial search and then the first 1000 points
	# output format: nsearch reduced-chiS Ndata stk dip rak dep lon lat t0
	is=1000; ie=999999
	Emul=1; chiS_max=999999
	awk 'BEGIN{flag=1;NFold=0}{if(NF==0&&NFold!=0){flag++} if(flag==1)print; NFold=NF}' $fMC | awk 'NF>0' | grep -v 'rejected' | awk -v is=$is -v ie=$ie '$2>is&&$2<ie' | awk -F\( '{print $3,$1}' | sed s/')'/''/ | awk -v Emul=$Emul '{print $15,$13*Emul,$10,$1,$2,$3,$4,$5,$6,$7}' | awk -v cm=$chiS_max '$2/$3<cm' > .PlotPosterior_tmp
	if [ `more .PlotPosterior_tmp | wc -l` -lt 5 ]; then
		echo "   Warning: problematic .PlotPosterior_tmp!"
		continue
	fi

	# location ellipse
	reslon=`ComputeMeanStd .PlotPosterior_tmp 8`
	reslat=`ComputeMeanStd .PlotPosterior_tmp 9`
	meanlon[idir]=`echo $reslon | awk '{lon=$1;if(lon<0.){lon+=360.}print lon}'`
	#lonavg=`echo $lonavg ${meanlon[idir]} | awk '{print $1+$2}'`
	stdlon[idir]=`echo $reslon | awk '{print $2}'`
	meanlat[idir]=`echo $reslat | awk '{print $1}'`
	stdlat[idir]=`echo $reslat | awk '{print $2}'`

	# best fitting model
	# check if the last line of fMC is marked as 'best'
	if [ `tail -n1 $fMC | grep -c 'best'` != 1 ]; then
		echo "best fitting model not found in "$fMC
		continue
	fi
	# tail -n1: best fitting from MC;  tail -n2: minimum from SA
	model_best[idir]=`tail -n1 .PlotPosterior_tmp | awk '{print $4,$5,$6,$7,$8,$9,$10}'`
	#cloc[idir]=`echo ${model_best[idir]} | awk '{print $5"_"$6}'`

	rm -f .PlotPosterior_tmp

	# extend region to plot
	ExtendRegion ${meanlon[idir]} ${meanlat[idir]} ${stdlon[idir]} ${stdlat[idir]}
	ExtendRegion `echo ${model_best[idir]} | awk '{print $5,$6}'`

	echo $rdir
	echo "   "${meanlon[idir]} ${meanlat[idir]}"   -   "${model_best[idir]}
	let idir++
done
#lonavg=`echo $lonavg $idir | awk '{print $1/$2}'`

### extend the region for locations from previous studies
for file in `ls previous_studies/results_*.txt`; do
	ExtendRegion `awk 'NR==1{print $5,$6}' $file`
done

#REG=`GetRegion cloc[@]`
lonmid=`echo $REG | awk -F/ '{print $1,$2}' | sed s/'-R'/''/ | awk '{print 0.5*($1+$2)}'`
echo $REG

### plot ###
gmtset HEADER_FONT_SIZE 15
gmtset LABEL_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE 10
gmtset PLOT_DEGREE_FORMAT ddd.xxxx
psout=$ename'_locations.ps'
pwd | psbasemap $REG -Jn${lonmid}/120. -Ba0.02f0.005WeSn -K -X4 -Y6 -P > $psout
PlotText $REG "Resulting Locations"
pscoast -R -J -N3/0/0/0 -W3 -S135/206/235 -O -K >> $psout

# plot the results (posterior as ellipses and best-fitting as stars)
ndir=$idir
labels=('R_1000' 'L_1000' 'B_1000' 'C_1000' 'R_500' 'L_500' 'B_500' 'C_500')
lnames=('Rayl' 'Love' 'Both' 'Both-Love' 'R_500' 'L_500' 'B_500' 'C_500')
Nlab=${#labels[@]}
#colorE=('lightred' 'forestgreen' 'steelblue' '60/60/60' 'green' 'red' 'blue' 'black')
#colorP=( '225/64/0' '65/100/128' ) #'60/60/60')
colors=( 'lightred' 'forestgreen' 'gold' 'steelblue' 'purple' 'black' )
for (( idir=0; idir<$ndir; idir++ )); do
	clocS=`echo ${model_best[idir]} | awk '{print $5,$6}'`
	label=`echo ${type[idir]} | cut -d_ -f1,2`
	# search for label name
	for (( ilab=0; ilab<$Nlab; ilab++ )); do
		if [ $label == ${labels[ilab]} ]; then	break; fi
	done
	if [ $ilab == $Nlab ]; then
		echo "Error: label "$label" not found in array \"labels\""
		exit
	fi
	# plot ellipse
	#PlotEllipse $psout ${meanlon[idir]} ${meanlat[idir]} ${stdlon[idir]} ${stdlat[idir]} $clocS ${colorE[ilab]}
	_color=${colors[idir]}
	PlotEllipse $psout ${meanlon[idir]} ${meanlat[idir]} ${stdlon[idir]} ${stdlat[idir]} $clocS $_color
	# point
	if [ `echo ${type[idir]} | cut -d_ -f3` == "Ei" ]; then
		_colorP=${colorP[0]}
	else
		_colorP=${colorP[1]} 
	fi
	#PlotPoint $psout $clocS $_colorP ${colorE[ilab]}
	#PlotPoint $psout $clocS $_color $_color
	PlotBeachball $psout `echo ${model_best[idir]} | awk '{print $1,$2,$3,$4,$5,$6}'` $_color $_color
done

# plot results from other studies
for file in `ls previous_studies/results_*.txt`; do
	study_name=`echo $file | awk -F/ '{print $NF}' | cut -d_ -f2 | cut -d. -f1`
	#PlotPoint $psout `awk 'NR==1{print $5,$6}' $file` black black
	PlotBeachball $psout `awk 'NR==1{print $1,$2,$3,$4,$5,$6}' $file` 80/80/80 80/80/80 $study_name
done

# legends
gmtset ANNOT_FONT_SIZE 12
dx1=0.15i
dx2=0.34i
gap=0.06i
size_a=0.3
size_e=0.6
pen=6
ulon=`echo $REG | sed s/'\-R'/''/ | awk -F/ '{print $2}'`
llat=`echo $REG | sed s/'\-R'/''/ | awk -F/ '{print $3}'`
echo $ulon $llat
pslegend -R -J -D${ulon}/${llat}/3./3.5/BR -G200 -F -O -K <<- EOF >> $psout
	G 0.03i
	S $dx1 e $size_e - $pen,${colors[0]} $dx2 `echo ${type[0]} | cut -d_ -f1,3 | sed s/'Ei'/'3D'/ | sed s/'B'/'Both'/ | sed s/'R'/'Rayl'/`
	G $gap
	S $dx1 e $size_e - $pen,${colors[1]} $dx2 `echo ${type[1]} | cut -d_ -f1,3 | sed s/'Ei'/'3D'/ | sed s/'B'/'Both'/ | sed s/'R'/'Rayl'/`
	G $gap
	S $dx1 e $size_e - $pen,${colors[2]} $dx2 `echo ${type[2]} | cut -d_ -f1,3 | sed s/'Ei'/'3D'/ | sed s/'B'/'Both'/ | sed s/'R'/'Rayl'/`
	G $gap
	S $dx1 e $size_e - $pen,${colors[3]} $dx2 `echo ${type[3]} | cut -d_ -f1,3 | sed s/'Ei'/'3D'/ | sed s/'B'/'Both'/ | sed s/'R'/'Rayl'/`
	G $gap
EOF
#legend -R -J -D${ulon}/${llat}/2.5/3.5/BR -G200 -F -O -K <<- EOF >> $psout
#	G 0.03i
#	S $dx1 a $size_a ${colorP[0]} $pen,${colorP[0]} $dx2 3D
#	G $gap
#	S $dx1 a $size_a ${colorP[1]} $pen,${colorP[1]} $dx2 1D
#	G 0.2i
#	S $dx1 e $size_e - $pen,${colorE[0]} $dx2 ${lnames[0]}
#	G $gap
#	S $dx1 e $size_e - $pen,${colorE[1]} $dx2 ${lnames[1]}
#	G $gap
#	S $dx1 e $size_e - $pen,${colorE[2]} $dx2 ${lnames[2]}
#	G $gap
#	S $dx1 e $size_e - $pen,${colorE[3]} $dx2 ${lnames[3]}

pwd | psxy -R -J -O >> $psout
echo $psout

