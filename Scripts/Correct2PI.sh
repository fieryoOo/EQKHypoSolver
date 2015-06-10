#!/bin/bash

if [ $# != 3 ]; then
	echo "Usage: "$0" [data-dis] [data-dir] [result-dir]"
	exit
fi

# input
dis=$1
datadir=$2
resdir=$3

# check directories
if [ ! -e $datadir ] || [ ! -e $resdir ]; then
	echo $datadir" or "$resdir" not found!"
	exit
fi

ftmp=fCorrect2PI.temp
for wtype in R L; do
	for fdata in `ls ${datadir}/${wtype}_Sta_grT_phT_Amp_*sec_dis${dis}.txt`; do
		per=`echo $fdata | awk -F/ '{print $NF}' | cut -d_ -f6 | sed s/'sec'/''/`
		fres=${resdir}/${wtype}_azi_data_pred_${per}sec.txt_sta
		echo $per $fdata $fres
		if [ ! -e $fres ]; then 
			echo "   Skipping: "$fres" not found."
			continue; 
		fi
		Nbeg=`awk 'BEGIN{N=0}{if(substr($1,1,1)=="#"){N=NR}}END{print N}' $fres`
		awk -v N=$Nbeg 'NR>N&&NF>0&&$10!=-12345' $fres | awk '{print $1,$2,$8}' > $ftmp
		/home/tianye/code/Programs/mini_tools/PickByLocation $ftmp $fdata ${ftmp}.out
		awk '{if(NF>7){print $4,$5,$6,$3,$8}else{print $2,$3,$4,$5,$6}}' ${ftmp}.out > ${ftmp}
		if [ `more ${ftmp} | wc -l` != `more $fdata | wc -l` ]; then 
			echo "   Skipping: failed!"
			continue; 
		fi
		mv $ftmp $fdata
	done
	rm -f ${ftmp} ${ftmp}.out
done
