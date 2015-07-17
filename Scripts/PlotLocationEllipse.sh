#!/bin/bash

PlotEllipse() {
	local _clon=$1
	local _clat=$2
}

### main ###
if [ $# != 2 ] && [ $# != 3 ] && [ $# != 4 ] && [ $# != 5 ]; then
   echo "Usage: "$0" [input_file] [starting_isearch] [ending_isearch (optional)] [chiS_max (optional)] [Emul (optional)]"
   exit
fi

### params
fin=$1
fps=${fin}.ps
is=$2
ie=1e10
if [ $# -ge 3 ]; then ie=$3; fi
chiS_max=999999
if [ $# -ge 4 ]; then chiS_max=$4; fi
Emul=1
if [ $# -ge 5 ]; then Emul=$5; fi


### discard results from the first search and grab qualified data from the second
ftmp=.PlotLocationEllipse_tmp
awk 'BEGIN{flag=1;NFold=0}{if(NF==0&&NFold!=0){flag++} if(flag==1)print; NFold=NF}' $fin | awk 'NF>0' | grep 'accepted' | awk -v is=$is -v ie=$ie '$2>is&&$2<ie' | awk -F\( '{print $3,$1}' | sed s/')'/''/ | awk -v Emul=$Emul '{print $6,$7,$14*Emul,$11,$16}' | awk -v cm=$chiS_max '$3/$4<cm' > $ftmp
#awk 'BEGIN{flag=1;NFold=0}{if(NF==0&&NFold!=0){flag++} if(flag==1)print; NFold=NF}' $fin | awk 'NF>0' | grep 'rejected' | awk -v is=$is -v ie=$ie '$2>is&&$2<ie' | awk -F\( '{print $3,$1}' | sed s/')'/''/ | awk -v Emul=$Emul '{print $16,$14*Emul,$11,$1,$2,$3,$4,$6,$7,$8}' > .PlotPosterior_tmp_rej


### check posterior file
Naccept=`more $ftmp | wc -l`
if [ $Naccept -le 0 ]; then
	echo "Problematic PosteriorD.txt"
	exit
fi
Nsearch=`tail -n1 $ftmp | awk '{if(NF>0){print $1}else{print 0}}'`


### exe for BivariateNormal fitting
exeF=/projects/yeti4009/code/Programs/mini_tools/BivariateNormal/FitBivariateNormal
if [ ! -e $exeF ]; then
	exeF=/home/tianye/code/Programs/mini_tools/BivariateNormal/FitBivariateNormal
fi
if [ ! -e $exeF ]; then
	echo "exe for BivariateNormal fitting not found!"
	exit
fi


### plotting starts

# plot settings
gmtset HEADER_FONT_SIZE 15
gmtset ANNOT_FONT_SIZE 10
gmtset HEADER_OFFSET 0.
gmtset PLOT_DEGREE_FORMAT ddd.xxxx

# region and psout
REG=`awk '{print $1,$2}' $ftmp | minmax -C | awk '{lonmin=$1;lonmax=$2;latmin=$3;latmax=$4;dlon=(lonmax-lonmin)*0.1;dlat=(latmax-latmin)*0.1;print lonmin-dlon, lonmax+dlon, latmin-dlat, latmax+dlat}' | awk '{print "-R"$1"/"$2"/"$3"/"$4}'`
clon=`awk '{print $1}' $ftmp | minmax -C | awk '{print 0.5*($1+$2)}'`
psout=LocationEllipse.ps

# plot locations
psxy $ftmp $REG -Jn${clon}/150 -Ba0.02f0.01:."Location distribution":WeSn -Sc0.03 -Gred -X4.5 -Y7 -K -P > $psout

# compute/plot error ellipse
prob=0.9
Einfo=`$exeF $ftmp ${prob} | awk '{print $1,$2,90.-$3,$4,$5}'`
text=`echo $Einfo | awk -v p=$prob '{printf "%.0f% : %.1fkm %.1fkm",p*100,$4*0.5,$5*0.5}'`
echo $Einfo
echo $Einfo | psxy -R -J -SE -W10,steelblue -O -K >> $psout

# ellipse info
loc=`awk '{print $1,$2}' $ftmp | minmax -C | awk '{lonmin=$1;lonmax=$2;latmin=$3;latmax=$4;dlon=(lonmax-lonmin)*0.1;dlat=(latmax-latmin)*0.1;print lonmin,latmax}'`
pstext -R -J -V -O -K -N << EOF >>  $psout
$loc 16 0.0 4 LT $text
EOF
#1.0 6 16 0.0 4 LT ${texts[1]}
#1.0 8 16 0.0 4 LT ${texts[2]}

pwd | psxy -R -J -O >> $psout
echo $psout

