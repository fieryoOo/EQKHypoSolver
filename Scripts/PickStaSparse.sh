#!/bin/bash

if [ $# != 3 ]; then
	echo "Usage: "$0" [clon] [clat] [#sta to pick]"
	exit
fi

### main ###
dis=1000
# list of fmeasurements to check
fcheck_list=()
fparam=param_base.txt
ifile=0
if [ -e $fparam ]; then
	for file in `grep -v "\#" param_base.txt  | awk '$1=="fRm"||$1=="fLm"' | awk '{print $2}' | sed s/'dis500'/'dis'$dis/g`; do
		fcheck_list[ifile]=$file
		let ifile++
	done
else
	echo -e "\n\n!!!Warning: "$fparam" not found, checking all measurement files!!!\n\n"
	for file in `ls Measurements/?_Sta_grT_phT_Amp_*sec_dis${dis}.txt`; do
		fcheck_list[ifile]=$file
		let ifile++
	done
fi

# check station list
fsta=../Scripts/station.lst
floc=station.loc; floctmp=temp.loc
awk '{print $2,$3,$1,0}' $fsta > $floc
for file in ${fcheck_list[@]}; do
	isgood=`/home/tianye/code/Programs/mini_tools/PickByLocation $floc $file $floctmp Y | awk 'NR>1{printf $1" "}END{print ""}' | awk '{if($1==$2){print "true"}else{print "false"}}'`
	if [ $isgood != true ]; then
		echo "station list "$fsta" not complete!"
		exit
	fi
done

# create list of valid stations
fstanew=./goodsta.lst
for file in ${fcheck_list[@]}; do
	echo $file
	/home/tianye/code/Programs/mini_tools/PickByLocation $file $floc $floctmp Y
	paste $floc $floctmp | awk '{N=$4; if(NF>6){N+=1} print $1,$2,$3,N}' > tempp
	mv tempp $floc
done

awk -v Nfile=${#fcheck_list[@]} '$4>=Nfile-1' $floc > tempp
mv tempp $floc
rm -f $floctmp

# compute azimuths
clon=$1
clat=$2
rm -f tempp
while read lon lat sta N; do
	azi=`/home/tianye/usr/bin/get_dist $clat $clon $lat $lon a1`
	echo $lon $lat $sta $azi >> tempp
done < $floc
sort -g -k4 tempp > $floc

# now randomly pick stations
exeRand=/home/tianye/usr/bin/Rand
Nsta=$3
binsize=`echo $Nsta | awk '{print 360./$1}'`
azil=`$exeRand 0 | awk -v binsize=$binsize '{print $1*binsize}'`
fresult=station_${clon}_${clat}_${Nsta}.txt
rm -f $fresult
for cazi in `echo $Nsta | awk -v azil=$azil -v Nsta=$Nsta '{for(i=0;i<Nsta;i++){cazi=azil+i*360./Nsta; print cazi}}'`; do
	rand=`$exeRand 1`
	azi=`echo $rand $cazi $binsize | awk '{print $2+$3*0.25*$1}'`
	awk -v azi=$azi '{if($4>azi){if(($4-azi)**2<(aziold-azi)**2){print $3,$1,$2}else{print staold,lonold,latold} exit}aziold=$4;lonold=$1;latold=$2;staold=$3}' $floc >> $fresult
done

echo $fresult
