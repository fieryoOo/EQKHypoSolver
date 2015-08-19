#!/bin/bash

if [ $# != 4 ] && [ $# != 6 ]; then
	echo "Usage: "$0" [clon] [clat] [#sta to pick] [picktype (0=even-azi, 1=even-space)] [max-dis (optional)] and [min-dis (optional)]"
	exit
fi

### main ###
clon=$1
clat=$2
Nsta=$3
picktype=$4
dis=1000; dismin=0
if [ $# == 6 ]; then
	dis=$5; dismin=$6
fi
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

# exe
#exePick=/home/tianye/code/Programs/mini_tools/PickByLocation
exePick=/projects/yeti4009/code/Programs/mini_tools/PickByLocation

# check station list
fsta=../Scripts/station.lst
floc=station.loc; floctmp=temp.loc
awk '{print $2,$3,$1,0}' $fsta > $floc
for file in ${fcheck_list[@]}; do
	isgood=`$exePick $floc $file $floctmp Y | awk 'NR>1{printf $1" "}END{print ""}' | awk '{if($1==$2){print "true"}else{print "false"}}'`
	if [ $isgood != true ]; then
		echo "station list "$fsta" not complete!"
		exit
	fi
done

# create list of valid stations
fstanew=./goodsta.lst
for file in ${fcheck_list[@]}; do
	awk -v dismin=$dismin '$6>dismin' $file > ${file}.tmp
	if [ `more ${file}.tmp | wc -l` == 0 ]; then
		echo "Warning: 0 valid station found in "$file
		echo "   Is distance info stored in column 6?   "
		continue
	fi
	$exePick ${file}.tmp $floc $floctmp Y 1> /dev/null
	rm -f ${file}.tmp
	paste $floc $floctmp | awk '{N=$4; if(NF>10){N+=1} print $1,$2,$3,N}' > tempp
	mv tempp $floc
done

awk -v Nfile=${#fcheck_list[@]} '$4>=Nfile-1' $floc > tempp
mv tempp $floc
rm -f $floctmp

# randomly pick stations from $floc
exeRand=/home/yeti4009/bin/Rand
exeDist=/home/yeti4009/bin/get_dist
fresult=station_${clon}_${clat}_${Nsta}.txt; rm -f $fresult
if [ $picktype == 0 ]; then
	# compute azimuths
	rm -f tempp
	while read lon lat sta N; do
		#azi=`$exeDist $clat $clon $lat $lon a1`
		azi=`$exeDist $clat $clon $lat $lon a1`
		echo $lon $lat $sta $azi >> tempp
	done < $floc
	sort -g -k4 tempp > temppp
	tail -n5 temppp | awk '{print $1,$2,$3,$4-360.}' > $floc
	cat temppp >> $floc
	head -n5 temppp | awk '{print $1,$2,$3,$4+360.}' >> $floc
	rm -f tempp temppp

	# now randomly pick stations
	#exeRand=/home/tianye/usr/bin/Rand
	binsize=`echo $Nsta | awk '{print 360./$1}'`
	azil=`$exeRand 0 | awk -v binsize=$binsize '{print $1*binsize}'`
	for cazi in `echo $Nsta | awk -v azil=$azil -v Nsta=$Nsta '{for(i=0;i<Nsta;i++){cazi=azil+i*360./Nsta; while(cazi<0.){cazi+=360.} while(cazi>=360.){cazi-=360.} print cazi}}'`; do
		rand=`$exeRand 1`
		azi=`echo $rand $cazi $binsize | awk '{print $2+$3*0.2*$1}'`
		line=`awk -v azi=$azi '{if($4>azi){if(($4-azi)**2<(aziold-azi)**2){print $0}else{print lineold} exit}aziold=$4;lineold=$0}' $floc`
		echo $line | awk '{print $3,$1,$2,$4}' >> $fresult
		awk -v line="$line" '$0!=line' $floc > tempp
		mv tempp $floc
		#awk -v azi=$azi '{if($4>azi){if(($4-azi)**2<(aziold-azi)**2){print $3,$1,$2,$4}else{print staold,lonold,latold,aziold} exit}aziold=$4;lonold=$1;latold=$2;staold=$3}' $floc >> $fresult
		#echo $azi $cazi $rand
	done
	#rm -f $floc
else
	isdis_min=70
	for ((i=0;i<$Nsta;i++)); do
		Nstaleft=`more $floc | wc -l`
		if [ $Nstaleft == 0 ]; then
			echo "Error: 0 station left!"
			exit
		fi
		# randomly pick
		ipick=`/home/yeti4009/bin/Rand 0 | awk -v N=$Nstaleft '{print int($1*N)+1}'`
		line=`awk -v i=$ipick 'NR==i{print $3,$1,$2,$4}' $floc` # $3=sta $1=lon $2=lat $4=azi
		echo $line >> $fresult
		# exclude nearby stations
		loc=`echo $line | awk '{print $3,$2}'`
		rm -f tempp
		while read lon lat sta azi; do
			isdis=`$exeDist $loc $lat $lon d`
			if [ `echo $isdis $isdis_min | awk '{if($1<$2){print 0}else{print 1}}'` == 1 ]; then
				echo $lon $lat $sta $azi >> tempp
			fi
		done < $floc
		mv tempp $floc
	done
fi

echo $fresult
