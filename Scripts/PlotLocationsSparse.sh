#!/bin/bash


PlotText() {
	local _REG=$1
	local _text=$2
	local _XS=$3
	local _YS=$4
	local _bbox=true
	if [ $# -ge 5 ]; then _bbox=$5; fi
	local _color=lightgray
	if [ $# -ge 6 ]; then _color=$6; fi

	local _flag=-W${_color},O3
	if [ $_bbox == false ]; then
		_flag="-G"${_color}
	fi

	#local _lb=('a' 'b' 'c' 'd' 'e' 'f')
	#local _title=`echo ${_text} | awk -v lb=${_lb[ifile]} '{print "("lb")  "$0}'`
	local _title=$_text
	echo ${_title} $_XS $_YS $_flag
	local _llon=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v xs=$_XS '{print $1+(0.01+xs)*($2-$1)}'`
	local _ulat=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v ys=$_YS '{print $4+(0.02+ys)*($4-$3)}'`
	echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "15. 0. 20 LT", $0}' | pstext -R -J $_flag -O -K -N >> $psout
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


ComputePlotEllipse() {
	local _fEdata=$1
	local _xtext=$2
	local _ytext=$3
	local _color=$4
	prob=0.9
	Einfo=`/projects/yeti4009/code/Programs/mini_tools/BivariateNormal/FitBivariateNormal $_fEdata $prob | awk '{print $1,$2,90.-$3,$4,$5}'`
   texts=`echo $Einfo | awk -v p=$prob '{printf "%.0f% : %.1fkm %.1fkm",p*100,$4*0.5,$5*0.5}'`
   echo $Einfo | psxy -R -J -SE -W5,${_color} -O -K >> $psout
	#echo "1.0 4 16 0.0 4 LT "$texts | pstext -R -J -V -X15 -Y0 -O -K -N >> $psout
	PlotText $REG "$texts" $_xtext $_ytext false $_color
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

	local _size1=0.2 #0.6
	local _size2=0.3 #0.65
	if [ `echo $_stk $_dip $_rak | awk '{if($1==0&&$2==0&&$3==0){print 0}else{print 1}}'` == 1 ]; then
		echo $_lon $_lat $_dep $_stk $_dip $_rak 5. 0 0 | psmeca -R -J -Sa${_size1} -G$_colorFill -W1,$_colorPen -O -K >> $_psout
		if [ $# -gt 8 ] && [ $_colorFill != $_colorPen ]; then
			echo $_lon $_lat | psxy -R -J -Sc${_size2} -W5,$_colorPen -O -K >> $_psout
		fi
	else
		echo $_lon $_lat | psxy -R -J -Sc${_size1} -G$_colorFill -W2,$_colorPen -O -K >> $_psout
	fi
	if [ $# -gt 9 ]; then
		echo $_lon $_lat $_text | awk '{print $1, $2-0.002, "8. 0. 20 LT", $3}' | pstext -R -J -Wlightgray,O3 -O -K -N >> $psout
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
fres=ResultsAll.txt; rm -f $fres
#for rdir in `ls -d results_SAMC_L_1000_1D_Sparse* results_SAMC_R_1000_1D_Sparse* results_SAMC_L_1000G_1D_Sparse* results_SAMC_B_1000_1D_Sparse* results_SAMC_L_1000_Ei_Sparse* results_SAMC_R_1000_Ei_Sparse* results_SAMC_L_1000G_Ei_Sparse* results_SAMC_B_1000_Ei_Sparse*`; do
for rdir in `ls -d results_SAMC_L_1000_1D_Sparse? results_SAMC_R_1000_1D_Sparse? results_SAMC_L_1000G_1D_Sparse? results_SAMC_B_1000_1D_Sparse? results_SAMC_L_1000_Ei_Sparse? results_SAMC_R_1000_Ei_Sparse? results_SAMC_L_1000G_Ei_Sparse? results_SAMC_B_1000_Ei_Sparse?`; do
#for rdir in `ls -d results_SAMC_B_1000_Ei_Sparse*`; do

	type[idir]=`echo $rdir | awk -F/ '{print $NF}' | cut -d_ -f3,4,5,6`

	### location range from MC (PosteriorD.txt)
	fMC=${rdir}/PosteriorD.txt
	if [ ! -e $fMC ]; then
		echo "   Warning: file "$fMC" not found"
		continue
	fi

if [ 0 == 1 ]; then
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
	ExtendRegion ${meanlon[idir]} ${meanlat[idir]} ${stdlon[idir]} ${stdlat[idir]}

	rm -f .PlotPosterior_tmp
fi

	# best fitting model
	# check if the last line of fMC is marked as 'best'
	line=`tail -n1 $fMC`
   if [ `echo $line | grep -c 'best'` != 1 ]; then
      echo "best fitting model not found in "$fMC
      continue
   fi
	echo $line | awk -F'minfo = \\(' '{print $2}' | awk -v type=${type[idir]} -F')' '{print $1,type}' >> $fres
	# this is wrong! model_best[idir]=`tail -n1 .PlotPosterior_tmp | awk '{print $4,$5,$6,$7,$8,$9,$10}'`
	model_best[idir]=`echo $line | awk -F'(' '{print $3}' | cut -d')' -f1 | awk '{print $1,$2,$3,$4,$5,$6,$7}'`
	#cloc[idir]=`echo ${model_best[idir]} | awk '{print $5"_"$6}'`


	# extend region to plot
	ExtendRegion `echo ${model_best[idir]} | awk '{print $5,$6}'`

	echo -e $rdir":\t"${model_best[idir]}
	#echo "   "${meanlon[idir]} ${meanlat[idir]}"   -   "${model_best[idir]}
	let idir++
done
#lonavg=`echo $lonavg $idir | awk '{print $1/$2}'`

### extend the region for locations from previous studies
predir=previous_studies
if [ ! -e $predir ]; then predir=../${predir}; fi
for file in `ls ${predir}/results_*.txt`; do
	ExtendRegion `awk 'NR==1{print $5,$6}' $file`
done

### input region? ###
if [ $# -gt 0 ];then
	REG=$1
fi

#REG=`GetRegion cloc[@]`
lonmid=`echo $REG | awk -F/ '{print $1,$2}' | sed s/'-R'/''/ | awk '{print 0.5*($1+$2)}'`
latmid=`echo $REG | awk -F/ '{print $3,$4}' | awk '{print 0.5*($1+$2)}'`
echo $REG

### plot ###
gmtset HEADER_FONT_SIZE 15
gmtset LABEL_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE 10
gmtset PLOT_DEGREE_FORMAT ddd.xxxx
psout=$ename'_locations.ps'
scale=`echo $REG | sed s/'-R'/''/ | awk -F/ '{ddeg1=18./($2-$1); ddeg2=20./($4-$3); if(ddeg1>ddeg2){ddeg1=ddeg2} printf "%.2g", ddeg1}'`
pwd | psbasemap $REG -Jn${lonmid}/${scale} -Ba0.1f0.02WeSn -K -X2.6 -Y3.5 -P > $psout
PlotText $REG "Resulting Locations"
pscoast -R -J -N3/0/0/0 -W3 -S135/206/235 -O -K >> $psout

# plot the results (posterior as ellipses and best-fitting as stars)
ndir=$idir
#labels=('R_1000' 'L_1000' 'B_1000' 'R_500' 'L_500' 'B_500')
#lnames=('Rayl' 'Love' 'Both' 'Both-Love' 'R_500' 'L_500' 'B_500' 'C_500')
labels=('B_1000_Ei' 'R_1000_Ei' 'L_1000_Ei' 'L_1000G_Ei' 'B_1000_1D' 'R_1000_1D' 'L_1000_1D' 'L_1000G_1D')
lnames=('3D_Both' '3D_Rayl' '3D_Love' '3D_LovG' '1D_Both' '1D_Rayl' '1D_Love' '1D_LovG')
colors=('black' 'brown' 'blue' 'forestgreen' 'darkgray' 'lightpink' 'lightblue' 'lightgreen')
Nlab=${#labels[@]}
#colorE=('lightred' 'forestgreen' 'steelblue' '60/60/60' 'red' 'green' 'blue' 'black')
#colorP=( '65/100/128' '225/64/0' ) #'60/60/60')
for (( idir=0; idir<$ndir; idir++ )); do
	clocS=`echo ${model_best[idir]} | awk '{print $5,$6}'`
	label=`echo ${type[idir]} | cut -d_ -f1,2,3`
	trlno=`echo ${type[idir]} | cut -d_ -f4 | sed s/'Sparse'/''/`
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
	_color=${colors[ilab]}
	#PlotEllipse $psout ${meanlon[idir]} ${meanlat[idir]} ${stdlon[idir]} ${stdlat[idir]} $clocS $_color
	# point
	#PlotPoint $psout $clocS $_colorP ${colorE[ilab]}
	#PlotPoint $psout $clocS $_color $_color
	PlotBeachball $psout `echo ${model_best[idir]} | awk '{print $1,$2,$3,$4,$5,$6}'` $_color #$_color $trlno
	PlotBeachball $psout `echo ${model_best[idir]} | awk '{print $1,$2,$3,$4,$5,$6}'` $_color #$_color $trlno
done

# plot results from other studies
for file in `ls ${predir}/results_*.txt`; do
	study_name=`echo $file | awk -F/ '{print $NF}' | cut -d_ -f2 | cut -d. -f1`
	#PlotPoint $psout `awk 'NR==1{print $5,$6}' $file` black black
	PlotBeachball $psout `awk 'NR==1{print $1,$2,$3,$4,$5,$6}' $file` orange orange $study_name
done


# compute/plot error ellipses
textloc=(0.05_0.12 0.05_0.09 0.05_0.06 0.05_0.03 0.55_0.12 0.55_0.09 0.55_0.06 0.55_0.03)
fEdata=EllipseData.tmp
ilab=0
for _label in ${labels[@]}; do
	if [ `grep $_label $fres | wc -l` == 0 ]; then continue; fi
	grep $_label $fres | awk '{print $5,$6,$8}' > $fEdata
	#echo "Check: "$_label $fEdata ${textloc[ilab]} ${colors[ilab]}
	ComputePlotEllipse $fEdata `echo ${textloc[ilab]} | awk -F_ '{print $1,$2}'` ${colors[ilab]}
	let ilab++
done
rm -f $fEdata


# legends
gmtset ANNOT_FONT_SIZE 12
dx1=0.15i
dx2=0.35i
gap=0.06i
size_a=0.25
size_e=0.6
pen=3
colpen=white
ulon=`echo $REG | sed s/'\-R'/''/ | awk -F/ '{print $2}'`
llat=`echo $REG | sed s/'\-R'/''/ | awk -F/ '{print $3}'`
#pslegend -R -J -D${ulon}/${llat}/3./3.5/BR -G200 -F -O -K <<- EOF >> $psout
#	G 0.03i
#	S $dx1 e $size_e - $pen,${colors[0]} $dx2 `echo ${type[0]} | cut -d_ -f1,3 | sed s/'Ei'/'3D'/ | sed s/'B'/'Both'/ | sed s/'R'/'Rayl'/`
#	G $gap
#	S $dx1 e $size_e - $pen,${colors[1]} $dx2 `echo ${type[1]} | cut -d_ -f1,3 | sed s/'Ei'/'3D'/ | sed s/'B'/'Both'/ | sed s/'R'/'Rayl'/`
#	G $gap
#	S $dx1 e $size_e - $pen,${colors[2]} $dx2 `echo ${type[2]} | cut -d_ -f1,3 | sed s/'Ei'/'3D'/ | sed s/'B'/'Both'/ | sed s/'R'/'Rayl'/`
#	G $gap
#	S $dx1 e $size_e - $pen,${colors[3]} $dx2 `echo ${type[3]} | cut -d_ -f1,3 | sed s/'Ei'/'3D'/ | sed s/'B'/'Both'/ | sed s/'R'/'Rayl'/`
#	G $gap
#EOF
echo $lonmid $latmid
pslegend -R -J -D${ulon}/${llat}/3./6.5/BR -G200 -F -O -K <<- EOF >> $psout
	G -0.15i
	M $lonmid $latmid 3 p
	G 0.12i
	S $dx1 c $size_a ${colors[0]} $pen,$colpen $dx2 ${lnames[0]}
	G $gap
	S $dx1 c $size_a ${colors[1]} $pen,$colpen $dx2 ${lnames[1]}
	G $gap
	S $dx1 c $size_a ${colors[2]} $pen,$colpen $dx2 ${lnames[2]}
	G $gap
	S $dx1 c $size_a ${colors[3]} $pen,$colpen $dx2 ${lnames[3]}
	G $gap
	S $dx1 c $size_a ${colors[4]} $pen,$colpen $dx2 ${lnames[4]}
	G $gap
	S $dx1 c $size_a ${colors[5]} $pen,$colpen $dx2 ${lnames[5]}
	G $gap
	S $dx1 c $size_a ${colors[6]} $pen,$colpen $dx2 ${lnames[6]}
	G $gap
	S $dx1 c $size_a ${colors[7]} $pen,$colpen $dx2 ${lnames[7]}
	G $gap
EOF

#echo 0.89 0.303 10 30 60 -110 5. 0 0 | psmeca -R0/1/0/1 -JX10 -Sa0.5 -G${colorP[0]} -W1,${colorP[0]} -N -O -K >> $psout
#echo 0.89 0.240 10 30 60 -110 5. 0 0 | psmeca -R0/1/0/1 -JX10 -Sa0.5 -G${colorP[1]} -W1,${colorP[1]} -N -O -K >> $psout

pwd | psxy -R -J -O >> $psout
echo $psout

