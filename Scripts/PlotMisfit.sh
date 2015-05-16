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
   #echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "15. 0. 20 LT", $0}' | pstext -R -J -Wlightgray,O3 -O -K -N >> $psout
   echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "12. 0. 20 LT", $0}' | pstext -R -J -O -K -N >> $psout
   let ifile++
}

GetSNR() {
   local _sta=$1
   local _f1=/work1/tianye/EQKLocation/SAC/Disps/${event}_disp_LHZ_10sec.txt
   local _f2=/work1/tianye/EQKLocation/SAC/Disps/${event}_disp_LHZ_16sec.txt

	#if [ $per == 10 ]; then
	#	SNR=`awk -v sta=$_sta '$4==sta{print $11}' $_f1`
	#else
	#	SNR=`awk -v sta=$_sta '$4==sta{print $11}' $_f2`
	#fi
   #SNR=`echo $SNR | awk '{ if(NF<1){ snr=2 }else{ snr=$1; if(snr<2){snr=2} } print snr} '`

   SNR1=`awk -v sta=$_sta '$4==sta{print $11}' $_f1`
   SNR2=`awk -v sta=$_sta '$4==sta{print $11}' $_f2`
   SNR=`echo $SNR1 $SNR2 | awk '{ if(NF<2){ snr=2 }else{ snr=$1; if(snr>$2){snr=$2} if(snr<2){snr=2} } print snr} '`
}

PlotHisto() {
	local _fin=$1
	local _REG=$2
	local _colorflag=$3
	local _ptext=$4
	local _YS=$5
	pshistogram $_fin $_REG -J ${_colorflag} -Z1 -O -K >> $psout
	if [ $_ptext != true ]; then return; fi
	local _avg=`awk 'BEGIN{a=0}{a+=$1}END{print a/NR}' $_fin`
	local _std=`awk -v avg=$_avg 'BEGIN{std=0}{std+=($1-avg)**2}END{print (std/(NR-1))**0.5}' $_fin`
	local _label=`echo $_avg $_std | awk '{printf "%.2f +/- %.2f\n",$1,$2}'`
	PlotText $_REG "$_label" 0.1 ${_YS}
}


### main ###
if [ $# != 1 ]; then
	echo "Usage: "$0" [result_file]"
	exit
fi

fin=$1
if [ ! -e $fin ]; then
	echo "bad file: "$fin
	exit
fi

# event
event=`pwd | awk -F/ '{print $(NF-1)}'`

# period
per=`echo $fin | sed s/'_'/' '/g | sed s/'\.'/' '/g | awk '{for(i=1;i<10;i++)print $i}' | grep sec | sed s/'sec'/''/`

# file label
#label=`readlink -f $fin | awk -F/ '{print $(NF-1)"_"$(NF)}' | cut -d. -f1 | awk -F_ '{print "data-used="$3" Map="$5" data-plot="$6" "$10}'`
label=`readlink -f $fin | awk -F/ '{print $(NF-1)"_"$(NF)}' | cut -d. -f1 | awk -F_ '{print $3" "$5" "$10}'`
echo $label

### set up gmt
gmtset HEADER_FONT_SIZE 12
gmtset LABEL_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE 10
gmtset HEADER_OFFSET 0.

### compute
Ntmp=`grep -n '#' $fin | awk -F: '{print $1}' | tail -n1`
#awk -v N=$Ntmp '{if(NR>N){print $1,$2,$8-$9-$10}}' $fin | awk 'NF==3' | psxy -R -J -Sc0.5 -W1,black -C$fcpt -O -K >> $psout
fdata=${fin}.tmp; rm -f $fdata
awk -v N=$Ntmp '{if(NR>N&&$9>-999&&$10>-999){print $1,$2,$8-$9,$10}}' $fin | awk 'NF==4' | while read lon lat Tmis Tsrc; do
#awk -v N=$Ntmp '{if(NR>N&&$9>-999&&$10>-999){print $1,$2,$5-$6,$7}}' $fin | awk 'NF==4' | while read lon lat Tmis Tsrc; do
	sta=`awk -v lon=$lon -v lat=$lat 'BEGIN{if(lon<0){lon+=360.}}{lonsta=$2; if(lonsta<0.){lonsta+=360.} latsta=$3; if( (lon-lonsta)**2+(lat-latsta)**2 < 0.0001 ){print $1} }' /home/tianye/EQKHypoSolver/Scripts/station.lst`
	GetSNR $sta
	#echo $sta $lon $lat $Tmis $SNR
	echo $lon $lat $Tmis $Tsrc $SNR | awk '{size=0.04*$5; if(size>0.6){size=0.6;} print $1,$2,$3,$4,size,$5}' >> $fdata
done

### plot
psout=${fin}.ps
REG=-R-123/-109/34/48

### with source
fcpt=/home/tianye/EQKHypoSolver/Scripts/Tsource.cpt
psbasemap $REG -Jn-116/0.6 -Ba5f1/a5f1:."Phase Time Misfit ($label)":WeSn -X4.5 -Y6 -K > $psout
pscoast -R -J -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 >> $psout
psxy /home/tianye/code/Programs/head/wus_province_II.dat -R -J -W5/255/0/0 -M"99999 99999"  -O -K >> $psout
psxy /home/tianye/code/Programs/head/platebound.gmt -R -J -W5/255/0/0 -M"99999 99999"  -O -K >> $psout
awk '{print $1,$2,$3,$5}' $fdata | psxy -R -J -Sc -W1,black -C$fcpt -O -K >> $psout
#psscale  -C$fcpt -B1:"Source Phase (sec)": -P -D4./-1./8./0.4h -O -K >> $psout
psscale  -C$fcpt -B:"Misfit (sec)": -P -D3.95/-1./7.9/0.4h -O -K >> $psout

# legend (SNR - size)
SNRs=( 3 8 15 )
sizes=( '' '' '' )
for ((i=0;i<3;i++)); do
   sizes[i]=`echo ${SNRs[i]} | awk '{size=0.04*$1; if(size>0.6){size=0.6;} print size}'`
done
echo ${sizes[@]} ${SNRs[@]}
pslegend -R -D-123.5/33.5/2.4/2.2/BL -J -F -G220 -O -K <<- EOF >> $psout
	G 0.05i
	S 0.1i c ${sizes[0]} white 0p 0.3i SNR=${SNRs[0]}
	G 0.1i
	S 0.1i c ${sizes[1]} white 0p 0.3i SNR=${SNRs[1]}
	G 0.1i
	S 0.1i c ${sizes[2]} white 0p 0.3i SNR=${SNRs[2]}
	G 0.05i
EOF


### final
fcpt=/home/tianye/EQKHypoSolver/Scripts/TshiftP.cpt
psbasemap -R -J -Ba5f1/a5f1:."Phase Time Misfit ($label)":WeSn -X9 -K -O >> $psout
pscoast -R -J -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 >> $psout
psxy /home/tianye/code/Programs/head/wus_province_II.dat -R -J -W5/255/0/0 -M"99999 99999"  -O -K >> $psout
psxy /home/tianye/code/Programs/head/platebound.gmt -R -J -W5/255/0/0 -M"99999 99999"  -O -K >> $psout
awk '{print $1,$2,$3-$4,$5}' $fdata | psxy -R -J -Sc -W1,black -C$fcpt -O -K >> $psout
psscale  -C$fcpt -B1.5:"Misfit (sec)": -P -D3.95/-1./7.9/0.4h -O -K >> $psout

pslegend -R -D-123.5/33.5/2.4/2.2/BL -J -F -G220 -O -K <<- EOF >> $psout
	G 0.05i
	S 0.1i c ${sizes[0]} white 0p 0.3i SNR=${SNRs[0]}
	G 0.1i
	S 0.1i c ${sizes[1]} white 0p 0.3i SNR=${SNRs[1]}
	G 0.1i
	S 0.1i c ${sizes[2]} white 0p 0.3i SNR=${SNRs[2]}
	G 0.05i
EOF


# plot histogram
REG=-R-8/8/0/40
ftmp_histo=ftmp_histo.txt
psbasemap $REG -JX6/8 -Ba2.f0.5/a10f2:."Phase Time Shift (sec)":WeSn -X9 -O -K >> $psout
awk '{if($6>8)print $3-$4}' $fdata > $ftmp_histo
PlotHisto $ftmp_histo $REG "-W0.4 -Gsteelblue -L3,steelblue3" true -0.05
awk '{if($6>8)print $3}' $fdata > $ftmp_histo
PlotHisto $ftmp_histo $REG "-W0.5 -L3,lightred" true -0.15

##awk '{if($6>10)print $3-$4}' $fdata | pshistogram -R-8/8/0/40 -JX6/8 -W0.4 -Gsteelblue -L3,steelblue3 -Z1 -X9 -O -K >> $psout
#awk '{if($6>10)print $3}' $fdata | pshistogram -R -J -W0.5 -Ba1.f0.2/a10f2:."Phase Time Shift (sec)":WeSn -L3,lightred -Z1 -O -K >> $psout
rm -f $fdata

### finalize
pwd | psxy -R -J -O >> $psout
echo $psout

