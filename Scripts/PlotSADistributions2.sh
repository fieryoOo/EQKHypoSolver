#!/bin/bash

PlotHisto() {
	local _fmod=$1
	local _columnN=$2
	local _XS=$3
	local _YS=$4
	local _label=$5
	local lb=$6
	local ub=$7
	local wd=$8
	local h_infile=.PlotPosterior_tmp
	is=0; ie=1e10; chiS_max=999999
	grep "search#" $_fmod | awk -v is=$is -v ie=$ie '$2>is&&$2<ie' | awk -F\( '{print $3,$1,$4}' | sed s/')'/''/g | awk -v Emul=$Emul -v is=$is '{if(NF==17){M0=$5}else{M0=1} printf "%d %f %d ",NR+is,$(NF-3)*Emul,$(NF-6); print $1,$2,$3,$4,M0,$(NF-11),$(NF-10),$(NF-9),$(NF);}' | awk -v cm=$chiS_max '$2/$3<cm' > ${h_infile}
	#grep 'rejected' .PlotPosterior_tmp > .PlotPosterior_tmp_rej
	#awk 'BEGIN{lastacc=""}{if($(NF)=="accepted"){print $0; lastacc=$0}else{if(lastacc!=""){print lastacc}}}' .PlotPosterior_tmp > .PlotPosterior_tmp_acc
   mean=`awk -v cN=$_columnN 'BEGIN{sum=0}{sum+=$cN}END{print sum/NR}' $h_infile`
   std=`awk -v cN=$_columnN -v mean=$mean 'BEGIN{sum=0}{sum+=($cN-mean)**2}END{std=sqrt(sum/(NR-1));if(std<0.0001)std=0.0001; print std}' $h_infile`
   std2=`echo $std | awk '{print $1*2}'`
	# region boundaries
	if [ $# == 5 ]; then
	   lb=`echo $mean $std | awk '{print $1-10.*$2}'`
		ub=`echo $mean $std | awk '{print $1+10.*$2}'`
	   wd=`echo $std | awk '{print $1*0.5}'`
	fi

	# plot histograms
   b_a=`echo $lb $ub | awk '{printf "%.1g",($2-$1)*0.2}'`
   b_f=`echo $b_a | awk '{print $1/5.}'`
   local REG=-R${lb}/${ub}/0/40
	_title=""
   awk -v cN=$_columnN '{print $cN}' $h_infile | pshistogram -JX4.5/5 $REG -Ba${b_a}f${b_f}-0.:"$_label":/a10f5:"No. Models (%)"::."$_title":WSn -W${wd} -L8,navyblue -Gsteelblue3 -Z1 -X${_XS} -Y${_YS}  -O -K >> ${fps}
   #awk -v cN=$_columnN '{print $cN}' $h_infile | pshistogram -J $REG -Ba${b_a}f${b_f}-0./a10f5:."$_title":WSn -W${wd} -Gsteelblue -L5,steelblue3 -Z1 -O -K >> ${fps}

	# mean and std
   #echo -e $mean 0"\n"$mean 100 | psxy -J -R -A -W8/255/100/100 -O -K >>${fps}
   l_2sig=`echo $mean $std2 | awk '{print $1-$2}'`
   u_2sig=`echo $mean $std2 | awk '{print $1+$2}'`
   echo -e $l_2sig 0"\n"$l_2sig 100 | psxy -J -R -A -W8/100/100/100 -O -K >>${fps}
   echo -e $u_2sig 0"\n"$u_2sig 100 | psxy -J -R -A -W8/100/100/100 -O -K >>${fps}
	echo $mean $std | awk '{mean=$1; std=$2; ftmp=log(2.0/std)/log(10); ord=ftmp==int(ftmp)?ftmp:(ftmp<0?int(ftmp):int(ftmp)+1); mul=10**ord; printf "0.1 9.9 10 0.0 4 LT %g (%g)",int(mean*mul+0.5)/mul,int(std*mul+0.5)/mul}' | pstext -R0/10/0/10 -Wwhite,o2 -J -V -O -K -N >> ${fps}
	#echo $mean $std | awk '{printf "0.1 9.9 10 0.0 4 LT %.7g +/- 2*%.2g",$1,$2}' | pstext -R0/10/0/10 -Wwhite,o2 -J -V -O -K -N >> ${fps}
	#echo -n $mean" +/- "$std"  "

	# plot reduced chi-square of the accepted models
	local REG=-R${lb}/${ub}/${rCl}/${rCu}
	local _ftmp=.PlotHisto_tmp
	awk -v cN=$_columnN '{print $cN, $2/$3}' $h_infile | sort -g -k1 > $_ftmp
	gmtset BASEMAP_FRAME_RGB 200/0/0
	psbasemap $REG -J -B/a0.5f0.1:"reduced chi-square":E -O -K >> ${fps}
	gmtset BASEMAP_FRAME_RGB +black
	#gmtset BASEMAP_FRAME_RGB black
	#psxy $_ftmp -R -J -A -Sc0.05 -Gred -O -K >> ${fps}
	psxy $_ftmp -R -J -A -W5,red -O -K >> ${fps}
	#local Nbin=2000; local Ntail=300
	#local Nb=100; local Nt=15
	#head -n $Ntail $_ftmp | awk 'BEGIN{min=99999}{if($2<min){x=$1;min=$2}}END{print x, min}' > ${_ftmp}_2
	#awk -v N=$Nbin 'BEGIN{min=99999}{if($2<min){x=$1;min=$2} if(NR%N==0){print x, min; min=99999}}' $_ftmp >> ${_ftmp}_2
	#tail -n $Ntail $_ftmp | awk 'BEGIN{min=99999}{if($2<min){x=$1;min=$2}}END{print x, min}' >> ${_ftmp}_2
	local Ntotal=`more $_ftmp | wc -l`
	local Ntail=`echo $Ntotal | awk '{print int($1/200+0.5)}'`; local Nt=`echo $Ntail | awk '{print int($1*0.05+0.5)}'`
	local Nbin=`echo $Ntotal $Ntail | awk '{print int(($1-2.*$2)*0.02)+1}'`; local Nb=`echo $Nbin | awk '{print int($1*0.02+0.5)}'`
	head -n $Ntail $_ftmp | sort -g -k2 | head -n $Nt > ${_ftmp}_2
	for ((ib=$Ntail+1;ib<$Ntotal-$Ntail;ib+=$Nbin)); do
		awk -v ib=$ib -v N=$Nbin 'NR>=ib&&NR<ib+N' $_ftmp | sort -g -k2 | head -n $Nb >> ${_ftmp}_2
	done
	tail -n $Ntail $_ftmp | sort -g -k2 | head -n $Nt >> ${_ftmp}_2
	#awk -v std=$std -v mean=$mean 'BEGIN{lb=mean-std*3.;ub=mean+std*3}{if($1>lb&&$1<ub)print}' ${_ftmp}_2 | psxy -R -J -A -Sc0.05 -O -K >> ${fps}
	#$exeBSpline ${_ftmp}_2 300.
	#awk -v std=$std -v mean=$mean 'BEGIN{lb=mean-std*3.;ub=mean+std*3}{if($1>lb&&$1<ub)print}' ${_ftmp}_2_spline | psxy -R -J -A -W6,220/0/0 -O -K >> ${fps}
	rm -f ${_ftmp} ${_ftmp}_2 ${_ftmp}_2_spline
}

###
#exeEx=./ExtractLowerModels2
exeEx=/lustre/janus_scratch/yeti4009/EQKLocation/TheoreticalResolution/ExtractLowerModels2
### main ###
if [ $# != 3 ]; then
	echo "Usage: "$0" [f_SAResults] [Emul] [rCu]"
	exit
fi

fdata=$1; Emul=$2;
rCl=0; rCu=$3;
frac=0.01;
fmod_stk=.fmod_stk.temp; fmod_dip=.fmod_dip.temp
fmod_rak=.fmod_rak.temp; fmod_dep=.fmod_dep.temp
if [ ! -e ${fmod_stk} ] || [ ! -e ${fmod_dip} ] || [ ! -e ${fmod_rak} ] || [ ! -e ${fmod_dep} ]; then
	${exeEx} $fdata $frac ${fmod_stk} ${fmod_dip} ${fmod_rak} ${fmod_dep}
	echo ${fmod_stk} ${fmod_dip} ${fmod_rak} ${fmod_dep}
fi

### plot histograms
gmtset HEADER_FONT_SIZE 12
gmtset LABEL_FONT_SIZE 10
gmtset ANNOT_FONT_SIZE 6
gmtset HEADER_OFFSET 0.
gmtset LABEL_OFFSET -0.15
gmtset ANNOT_OFFSET 0.05


fps=${fdata}.ps
pwd | psxy -R0/10/0/10 -JX1 -K -P > ${fps}
#pstext -R0/10/0/10 -JX15/25.5 -Wwhite,O3 -V -O -K -N << EOF >>  ${fps}
#2.8 10. 16 0.0 4 LT Naccept/Nsearch = $Naccept/$Nsearch
#EOF

PlotHisto ${fmod_stk} 4 2 16 "strike (deg)" 0 360 1
PlotHisto ${fmod_dip} 5 7 0 "dip (deg)" 0 90 1
PlotHisto ${fmod_rak} 6 -7 -6.5 "rake (deg)" -180 180 1
PlotHisto ${fmod_dep} 7 7 0 "depth (km)" 0 60 1

pwd | psxy -R -J -O >> ${fps}
rm -f .PlotPosterior_tmp .PlotPosterior_tmp_acc .PlotPosterior_tmp_rej
#rm -f ${fmod_stk} ${fmod_dip} ${fmod_rak} ${fmod_dep}

echo ${fps}

