#!/bin/bash

Alpha() {
	local per_v=$1
	local wtype=$2
	local alpha=-12345.
	if [ $wtype == 'R' ]; then
		if [ $per_v -ge 7 ] && [ $per_v -lt 9 ]; then
			v=2.6; Q=80	
		elif [ $per_v -ge 9 ] && [ $per_v -lt 11 ]; then
		   v=2.8; Q=80 # 102 (need to be decreased?)
		elif [ $per_v -ge 11 ] && [ $per_v -lt 13 ]; then
		   v=2.83; Q=81 
		elif [ $per_v -ge 13 ] && [ $per_v -lt 15 ]; then
		   v=2.86; Q=83
		elif [ $per_v -ge 15 ] && [ $per_v -lt 17 ]; then
		   v=2.9; Q=85
		elif [ $per_v -ge 17 ] && [ $per_v -lt 20 ]; then
		   v=2.95; Q=87
		elif [ $per_v -ge 20 ] && [ $per_v -lt 25 ]; then
		   v=3.0; Q=90
		elif [ $per_v -ge 25 ] && [ $per_v -lt 35 ]; then
		   v=3.25; Q=100
		elif [ $per_v -ge 35 ] && [ $per_v -lt 45 ]; then
		   v=3.5; Q=120	# 112 (need to be increased?)
	   else
		   echo "Error(AlphaR): alpha undefined for per = "$per_v"!"
			exit
	   fi
	elif [ $wtype == 'L' ]; then
		if [ $per_v -ge 7 ] && [ $per_v -lt 9 ]; then
			v=2.8; Q=80
		elif [ $per_v -ge 9 ] && [ $per_v -lt 11 ]; then
		   v=3.1; Q=80
		elif [ $per_v -ge 11 ] && [ $per_v -lt 13 ]; then
		   v=3.15; Q=81
		elif [ $per_v -ge 13 ] && [ $per_v -lt 15 ]; then
		   v=3.21; Q=83
		elif [ $per_v -ge 15 ] && [ $per_v -lt 17 ]; then
		   v=3.27; Q=85
		elif [ $per_v -ge 17 ] && [ $per_v -lt 20 ]; then
		   v=3.32; Q=87
		elif [ $per_v -ge 20 ] && [ $per_v -lt 25 ]; then
		   v=3.38; Q=90
		elif [ $per_v -ge 25 ] && [ $per_v -lt 35 ]; then
		   v=3.53; Q=100
		elif [ $per_v -ge 35 ] && [ $per_v -lt 45 ]; then
		   v=3.73; Q=120
	   else
		   echo "Error(AlphaR): alpha undefined for per = "$per_v"!"
			exit
	   fi
	else
		echo "Error(Alpha): unknown surface wave type: "$wtype"!!"
		exit
	fi
	alpha=`echo $per_v $v $Q | awk '{pi=3.1415926536; print pi/($1*$2*$3)}'`
	echo $alpha $Q
}


### main ###
if [ $# != 3 ]; then
	echo "Usage: "$0" [event name] [max dist] [snr min]"
	exit
fi

# parameters
event=`echo $1 | awk -F/ '{print $1}'`
dismax=$2
snrmin=$3
if [ $dismax -le 0 ] || [ $dismax -gt 999999 ]; then
	echo "Error(main): Invalid dis max: "$dismax
	exit
fi
dir_disp='SAC/Disps'
dir_out=${event}/Measurements

mkdir -p ${event}
mkdir -p ${dir_out}

### Extract R/L wave data from Z/T component ###
fQ=${dir_out}/Qvalues.txt
rm -f $fQ
correctGS=false
correctQ=false
echo "dismax = "$dismax >> $fQ
echo "snrmin = "$snrmin >> $fQ
echo "correctGS? : "$correctGS >> $fQ
echo "correctQ ? : "$correctQ  >> $fQ
for wtype in R L; do
	# component (file) name
	if [ $wtype == 'R' ]; then
		cpnt='Z'
	elif [ $wtype == 'L' ]; then
		cpnt='T'
	else
		echo "Error(main): unknown SW type: "$wtype
		exit
	fi
	# iterate through each file
	for file in `ls ${dir_disp}'/'${event}'_disp_LH'${cpnt}'_'*'sec.txt'`; do
		per=`echo $file | cut -d_ -f4 | cut -d. -f1`
	   per_v=`echo $per | sed s/'sec'/''/`
		# get attenuation coefficient
		alp_Q=`Alpha $per_v ${wtype}`
		echo $wtype $per_v ${alp_Q} >> $fQ
		if [ $correctQ == true ]; then
			alpha=`echo ${alp_Q} | awk '{print $1}'`
		else
			alpha=0
		fi
		fout=${dir_out}/${wtype}_Sta_grT_phT_Amp_${per}_dis${dismax}.txt
		# get lon-lat-grv-phv-amp info (amp is corrected for both geometric spreading and attenuation)
	   awk -v cGS=$correctGS -v per=$per_v -v alp=$alpha -v dismax=$dismax -v snrmin=$snrmin 'NR>1&&$11>snrmin&&$12>0&&$13>0&&$14>0&&$7<dismax{
			if( cGS=="false" ) {
				disQ = 1000.;
			} else {
				disQ = $7;
			}
			print $5,$6, $13, $14, $12*sqrt(disQ/1000.)/exp(-alp*$7)
		}' $file > $fout
		echo $file $fout
	done # file
done # wtype
