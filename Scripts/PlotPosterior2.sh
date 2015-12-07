#!/bin/bash

PlotHisto() {
   h_infile=$1
   columnN=$2
   Xoffset=$3
   Yoffset=$4
	local _title=$5
	local _label=$6
   mean=`awk -v cN=$columnN 'BEGIN{sum=0}{sum+=$cN}END{print sum/NR}' $h_infile`
   std=`awk -v cN=$columnN -v mean=$mean 'BEGIN{sum=0}{sum+=($cN-mean)**2}END{print sqrt(sum/(NR-1))}' $h_infile`
   std2=`echo $std | awk '{print $1*2}'`
	# region boundaries
   lb=`echo $mean $std | awk '{print $1-5.*$2}'`
   ub=`echo $mean $std | awk '{print $1+5.*$2}'`
   wd=`echo $std | awk '{print $1*0.8}'`

	# plot histograms
   b_a=`echo $std2 | awk '{printf "%.1g",$1}'`
   b_f=`echo $b_a | awk '{print $1/5.}'`
   local REG=-R${lb}/${ub}/0/40
	_title=""
   awk -v cN=$columnN '{print $cN}' $h_infile | pshistogram -JX4.5/5 $REG -Ba${b_a}f${b_f}-0.:"$_label":/a10f5:"No. Models (%)"::."$_title":WSn -W${wd} -L8,navyblue -Gsteelblue3 -Z1 -X${Xoffset} -Y${Yoffset}  -O -K >> ${fps}
   #awk -v cN=$columnN '{print $cN}' $h_infile | pshistogram -J $REG -Ba${b_a}f${b_f}-0./a10f5:."$_title":WSn -W${wd} -Gsteelblue -L5,steelblue3 -Z1 -O -K >> ${fps}

	# mean and std
   #echo -e $mean 0"\n"$mean 100 | psxy -J -R -A -W8/255/100/100 -O -K >>${fps}
   l_2sig=`echo $mean $std2 | awk '{print $1-$2}'`
   u_2sig=`echo $mean $std2 | awk '{print $1+$2}'`
   echo -e $l_2sig 0"\n"$l_2sig 100 | psxy -J -R -A -W8/100/100/100 -O -K >>${fps}
   echo -e $u_2sig 0"\n"$u_2sig 100 | psxy -J -R -A -W8/100/100/100 -O -K >>${fps}
	echo $mean $std | awk '{mean=$1; std=$2; ftmp=log(2.0/std)/log(10); ord=ftmp==int(ftmp)?ftmp:(ftmp<0?int(ftmp):int(ftmp)+1); mul=10**ord; printf "0.1 9.9 10 0.0 4 LT %g (%g)",int(mean*mul+0.5)/mul,int(std*mul+0.5)/mul}' | pstext -R0/10/0/10 -Wwhite,o2 -J -V -O -K -N >> ${fps}
	#echo $mean $std | awk '{printf "0.1 9.9 10 0.0 4 LT %.7g +/- 2*%.2g",$1,$2}' | pstext -R0/10/0/10 -Wwhite,o2 -J -V -O -K -N >> ${fps}
	echo -n $mean" +/- "$std"  "

	# plot reduced chi-square of the accepted models
	local REG=-R${lb}/${ub}/${rCl}/${rCu}
	local _ftmp=.PlotHisto_tmp
	awk -v cN=$columnN '{print $cN, $2/$3}' $h_infile | sort -g -k1 > $_ftmp
	gmtset BASEMAP_FRAME_RGB 200/0/0
	psbasemap $REG -J -B/a0.1f0.02:"reduced chi-square":E -O -K >> ${fps}
	gmtset BASEMAP_FRAME_RGB +black
	#gmtset BASEMAP_FRAME_RGB black
	#psxy $_ftmp -R -J -A -Sc0.01 -Gred -O -K >> ${fps}
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
	$exeBSpline ${_ftmp}_2 300.
	#$exeCSpline ${_ftmp}_2
	#awk -v std=$std -v mean=$mean 'BEGIN{lb=mean-std*3.;ub=mean+std*3}{if($1>lb&&$1<ub)print}' ${_ftmp}_2 | psxy -R -J -A -Sc0.05 -O -K >> ${fps}
	awk -v std=$std -v mean=$mean 'BEGIN{lb=mean-std*3.;ub=mean+std*3}{if($1>lb&&$1<ub)print}' ${_ftmp}_2_spline | psxy -R -J -A -W6,220/0/0 -O -K >> ${fps}
	rm -f ${_ftmp} ${_ftmp}_2 ${_ftmp}_2_spline

}

### main ###
if [ $# != 3 ] && [ $# != 4 ] && [ $# != 5 ]; then
   echo "Usage: "$0" [input_file] [Emul] [starting_isearch] [ending_isearch (optional)] [chiS_max (optional)]"
   exit
fi

### params
fin=$1
fps=${fin}.ps
Emul=$2
is=$3
ie=1e10
if [ $# -ge 4 ]; then ie=$4; fi
chiS_max=999999
if [ $# -ge 5 ]; then chiS_max=$5; fi
exeBSpline=/projects/yeti4009/code/Programs/Splines/BSpline1D
if [ ! -e $exeBSpline ]; then
	exeBSpline=/home/tianye/code/Programs/Splines/BSpline1D
fi
exeCSpline=/projects/yeti4009/code/Programs/Splines/CubicSpline1D
if [ ! -e $exeCSpline ]; then
	exeCSpline=/home/tianye/code/Programs/Splines/CubicSpline1D
fi

### discard results from the first search and grab qualified data from the second
awk 'BEGIN{flag=1;NFold=0}{if(NF==0&&NFold!=0){flag++} if(flag==1)print; NFold=NF}' $fin | awk 'NF>0' | awk -v is=$is -v ie=$ie '$2>is&&$2<ie' | awk -F\( '{print $3,$1,$4}' | sed s/')'/''/g | awk -v Emul=$Emul '{if(NF==17){M0=$5}else{M0=1} printf "%d %f %d ",$(NF-1),$(NF-3)*Emul,$(NF-6); print $1,$2,$3,$4,M0,$(NF-11),$(NF-10),$(NF-9),$(NF);}' | awk -v cm=$chiS_max '$2/$3<cm' > .PlotPosterior_tmp
grep 'accepted' .PlotPosterior_tmp > .PlotPosterior_tmp_acc
grep 'rejected' .PlotPosterior_tmp > .PlotPosterior_tmp_rej
#awk 'BEGIN{flag=1;NFold=0}{if(NF==0&&NFold!=0){flag++} if(flag==1)print; NFold=NF}' $fin | awk 'NF>0' | grep 'rejected' | awk -v is=$is -v ie=$ie '$2>is&&$2<ie' | awk -F\( '{print $3,$1}' | sed s/')'/''/ | awk -v Emul=$Emul '{print $16,$14*Emul,$11,$1,$2,$3,$4,$6,$7,$8}' > .PlotPosterior_tmp_rej
#awk 'NF!=0' $fin | awk 'begin{flag=0;Nold=0}{if($1<Nold){flag=1} Nold=$1; if(flag==1)print}' | awk -v is=$is -v ie=$ie '$1>is&&$1<ie' > .PlotPosterior_tmp


### check posterior file
Naccept=`more .PlotPosterior_tmp_acc | wc -l`
if [ $Naccept -le 0 ]; then
	echo "Problematic PosteriorD.txt"
	exit
fi
#Nsearch=`tail -n1 .PlotPosterior_tmp_acc | awk '{if(NF>0){print $1}else{print 0}}'`
Nsearch=`tail -n3 $fin | awk 'BEGIN{N=0}{if($2>N){N=$2}}END{print N}'`

### range for reduced chi-square
rCrange=`awk '{print $2/$3}' .PlotPosterior_tmp_acc | minmax -C`
rCl=`echo $rCrange | awk '{printf "%.2f", $1-0.02}'`
rCu=`echo $rCl | awk '{printf "%.2f", $1+0.2}'`
#rCl=2.1; rCu=2.4

### plot histograms
gmtset HEADER_FONT_SIZE 12
gmtset LABEL_FONT_SIZE 10
gmtset ANNOT_FONT_SIZE 6
gmtset HEADER_OFFSET 0.
gmtset LABEL_OFFSET -0.15
gmtset ANNOT_OFFSET 0.05

pwd | psxy -R0/1/0/1 -JX1 -K -P > ${fps}
pstext -R0/10/0/10 -JX15/25.5 -Wwhite,O3 -V -O -K -N << EOF >>  ${fps}
2.8 9.2 16 0.0 4 LT Naccept/Nsearch = $Naccept/$Nsearch
EOF
#4.0 10.2 16 0.0 4 LT Naccept = $Naccept

rm -f results

PlotHisto .PlotPosterior_tmp_acc 9 -0.9 17 "longitude" "longitude (deg)" >> results
PlotHisto .PlotPosterior_tmp_acc 10 6.7 0 "latitude" "latitude (deg)" >> results
PlotHisto .PlotPosterior_tmp_acc 11 6.7 0 "origin time" "origin time (sec)" >> results

PlotHisto .PlotPosterior_tmp_acc 4 -13.4 -6.5 "strike" "strike (deg)" >> results
PlotHisto .PlotPosterior_tmp_acc 5 6.7 0 "dip" "dip (deg)" >> results
PlotHisto .PlotPosterior_tmp_acc 6 6.7 0 "rake" "rake (deg)" >> results
PlotHisto .PlotPosterior_tmp_acc 7 -10.7 -6.2 "depth" "depth (km)" >> results
PlotHisto .PlotPosterior_tmp_acc 8 8 0 "M0" "M0" >> results

if [ 0 == 1 ]; then

finfo1=`more results | awk '{print $1,$4,$7,$10}'`
finfo2=`/projects/yeti4009/eqkhyposolver/Auxiliary $finfo1`
#finfo2=`/home/tianye/EQKHypoSolver/Auxiliary $finfo1`
echo -e "\nFocal: ("$finfo1") - ("$finfo2")" >> results

# plot misfits
# rejected
ie=`awk 'BEGIN{imax=0}{if(imax<$1){imax=$1}}END{print imax}' .PlotPosterior_tmp_acc`
awk '{print $1,$2/$3}' .PlotPosterior_tmp_rej | psxy -R${is}/${ie}/${rCl}/${rCu} -JX18/5 -Ba10000f2000/a0.1f0.02:."reduced chi-square":WeSn -A -Sc0.03 -G100/100/100 -X-13.7 -Y-6.5 -O -K >> ${fps}
# accepted
awk '{print $1,$2/$3}' .PlotPosterior_tmp_acc | psxy -R -J -A -Sc0.03 -Gred -O -K >> ${fps}
# accepted min boundary (min each 1000 acceptances)
awk '{print $1,$2/$3}' .PlotPosterior_tmp_acc |  awk 'BEGIN{min=99999;N=0;issum=0}{if(N==1000 && min!=99999){print issum/N,min; min=99999;issum=0;N=0}if($2<min){min=$2} N++;issum+=$1}' | psxy -R -J -A -S-1. -W5,red -O -K >> ${fps}

fi


pwd | psxy -R -J -O >> ${fps}
rm -f .PlotPosterior_tmp .PlotPosterior_tmp_acc .PlotPosterior_tmp_rej

echo ${fps}
