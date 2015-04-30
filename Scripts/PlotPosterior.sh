#!/bin/bash

PlotHisto() {
   h_infile=$1
   columnN=$2
   Xoffset=$3
   Yoffset=$4
	local _title=$5
   mean=`awk -v cN=$columnN 'BEGIN{sum=0}{sum+=$cN}END{print sum/NR}' $h_infile`
   std=`awk -v cN=$columnN -v mean=$mean 'BEGIN{sum=0}{sum+=($cN-mean)**2}END{print sqrt(sum/(NR-1))}' $h_infile`
   std2=`echo $std | awk '{print $1*2}'`
   lb=`echo $mean $std | awk '{print $1-5.*$2}'`
   ub=`echo $mean $std | awk '{print $1+5.*$2}'`
   wd=`echo $std | awk '{print $1*0.8}'`
   b_a=`echo $std2 | awk '{printf "%.1g",$1}'`
   b_f=`echo $b_a | awk '{print $1/5.}'`
   REG=-R${lb}/${ub}/0/40
   awk -v cN=$columnN '{print $cN}' $h_infile | pshistogram -JX5/6 $REG -Ba${b_a}f${b_f}-0./a10f5:."$_title":WeSn -W${wd} -Gsteelblue -L5,steelblue3 -Z1 -X${Xoffset} -Y${Yoffset} -O -K >> ${fps}
   echo -e $mean 0"\n"$mean 100 | psxy -J -R -A -W8/255/100/100 -O -K >>${fps}
   l_2sig=`echo $mean $std2 | awk '{print $1-$2}'`
   u_2sig=`echo $mean $std2 | awk '{print $1+$2}'`
   echo -e $l_2sig 0"\n"$l_2sig 100 | psxy -J -R -A -W8/100/100/100 -O -K >>${fps}
   echo -e $u_2sig 0"\n"$u_2sig 100 | psxy -J -R -A -W8/100/100/100 -O -K >>${fps}
	echo $mean $std | awk '{printf "1. 9.5 10 0.0 4 LT %.7g +/- 2*%.2g",$1,$2}' | pstext -R0/10/0/10 -J -V -O -K -N >> ${fps}
	echo -n $mean" +/- "$std"  "
}

### main ###
if [ $# != 2 ] && [ $# != 3 ]; then
   echo "Usage: "$0" [input_file] [starting_isearch] [ending_isearch (optional)]"
   exit
fi
fin=$1
fps=${fin}.ps
is=$2
ie=1e10
if [ $# == 3 ]; then ie=$3; fi
### discard results from the first search and grab qualified data from the second
awk 'NF>0' $fin | awk 'BEGIN{flag=1;Nold=0}{if($2<Nold){flag++} if(flag==2)print; Nold=$2}' | grep 'accepted' $fin | awk -v is=$is -v ie=$ie '$2>is&&$2<ie' | awk -F\( '{print $3,$1}' | sed s/')'/''/ | awk '{print $15,$1,$2,$3,$4,$5,$6,$7}' > .PlotPosterior_tmp
#awk 'NF!=0' $fin | awk 'begin{flag=0;Nold=0}{if($1<Nold){flag=1} Nold=$1; if(flag==1)print}' | awk -v is=$is -v ie=$ie '$1>is&&$1<ie' > .PlotPosterior_tmp
Naccept=`more .PlotPosterior_tmp | wc -l`
if [ $Naccept -le 0 ]; then
	echo "Problematic PosteriorD.txt"
	exit
fi
Nsearch=`tail -n1 .PlotPosterior_tmp | awk '{if(NF>0){print $1}else{print 0}}'`

gmtset HEADER_OFFSET -0.6
gmtset ANNOT_FONT_SIZE 6
pwd | psxy -R0/1/0/1 -JX1 -K -P > ${fps}
pstext -R0/10/0/10 -JX15/25.5 -V -O -K -N << EOF >>  ${fps}
1.0 9.8 16 0.0 4 LT Nsearch = $Nsearch
4.0 9.8 16 0.0 4 LT Naccept = $Naccept
EOF

rm -f results
PlotHisto .PlotPosterior_tmp 2 0 17.5 "strike" >> results
PlotHisto .PlotPosterior_tmp 3 8 0 "dip" >> results
PlotHisto .PlotPosterior_tmp 4 -8 -7.5 "rake" >> results
PlotHisto .PlotPosterior_tmp 5 8 0 "depth" >> results

PlotHisto .PlotPosterior_tmp 6 -9.5 -8.5 "longitude" >> results
PlotHisto .PlotPosterior_tmp 7 7 0 "latitude" >> results
PlotHisto .PlotPosterior_tmp 8 7 0 "origin time" >> results

finfo1=`more results | awk '{print $1,$4,$7,$10}'`
#finfo2=`/projects/yeti4009/code/Programs/ExploitEvent/SearchLocation/Auxiliary $finfo1`
finfo2=`/home/tianye/EQKHypoSolver/Auxiliary $finfo1`
echo -e "\nFocal: ("$finfo1") - ("$finfo2")" >> results

pwd | psxy -R -J -O >> ${fps}
rm -f .PlotPosterior_tmp

echo ${fps}
