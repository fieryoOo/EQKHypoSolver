#!/bin/bash

RewriteParams() {
	local _fin=$1
	local _fou=$2
	echo $refloc | awk '{printf "lon %f\nlat %f\n", $2, $1}' > $_fou
	echo $reffoc1 | awk '{printf "stk %f\ndip %f\nrak %f\ndep %f\n", $1,$2,$3,$4}' >> $_fou
	awk '$1!="lon" && $1!="lat" && $1!="stk" && $1!="dip" && $1!="rak" && $1!="dep"' $_fin >> $_fou
}

if [ $# != 1 ]; then
	echo "Usage: "$0" [res dir] OR [Rerunning.list for moving res directories back]"
	exit
fi

resdir=$1
if [ -f $resdir ]; then
	# Input is Rerunning.list!
	allsuc=true;
	while read dirnew dirold stmp; do
		if [ ! -e $dirnew ] || [ -e $dirold ]; then
			echo "Warning: new res "$dirnew" not exists / old res dir "$dirold" exists!"
			allsuc=false; continue
		fi
		echo $dirnew" -> "$dirold
		mv $dirnew $dirold
	done < $resdir
	if [ allsuc == true ]; then rm -f $resdir; fi
	exit
elif [ ! -d $resdir ]; then
	echo $resdir" is not a directory nor a regular file!"
	exit
fi

#minfo = (266.747  36.490  -42.753  5.952 2.216e+24   245.1035 41.1460  0.8630)	N = 15	E = 1.4747 (best)
refloc="41.145 245.09"
reffoc1="256. 36. -48. 5.5"
reffoc2="28. 64. -116. 5.5"
uncfoc="7. 3. 6. 0.5"

list_rerun=Rerunning.list
for dir in `ls -d $resdir/results_SAMC_B_*Ei* $resdir/results_SAMC_R_*Ei*`; do
	# find best fitting model
	fPos=${dir}/PosteriorD.txt
	if [ ! -e $fPos ]; then
		echo $fPos not found!!!
		continue
	fi
	line=`tail -n1 $fPos`
	if [ `echo $line | grep -c best` == 0 ]; then
		echo $fPos" incomplete: best fitting model not found at the end!!"
		continue
	fi
	model_best=`echo $line | awk -F"minfo = " '{print $2}' | awk -F"N = " '{print $1}' | sed s/'('/''/ | sed s/')'/''/`

	# analyse result
	E=`echo $line | awk -F"E = " '{print $2}' | awk '{print $1}'`
	loc=`echo $model_best | awk '{print $(NF-1),$(NF-2)}'`
	foc=`echo $model_best | awk '{print $1,$2,$3,$4}'`
	misloc=`/home/yeti4009/bin/get_dist $refloc $loc d`
	misfoc1=`echo $foc $reffoc1 $uncfoc | awk 'BEGIN{mis=0}{for(i=1;i<5;i++){ misC=(($(i)-$(i+4))/$(i+8));if(misC<0){misC=-misC} mis+=misC}}END{print mis}'`
	misfoc2=`echo $foc $reffoc2 $uncfoc | awk 'BEGIN{mis=0}{for(i=1;i<5;i++){ misC=(($(i)-$(i+4))/$(i+8));if(misC<0){misC=-misC} mis+=misC}}END{print mis}'`
	misfoc=`echo $misfoc1 $misfoc2 | awk '{if($1<$2){print $1}else{print $2}}'`
	if [ `echo $misfoc | awk '{if($1<10.){print 1}else{print 0}}'` == 1 ]; then continue; fi
	#echo $fPos $model_best "E = "$E" misloc = "$misloc" misfoc = "$misfoc

	# result is bad: relocate/rerun
	dirnew=`echo $dir | awk -F/ '{print "./"$(NF)}'`
	if [ -e $dirnew ]; then
		echo "targeting directory "$dirnew" exists!!!"
		exit
	fi
	if [ -e ${dir}"~" ]; then
		echo "saving directory "${dir}"~ exists!!!"
		exit
	fi
	mkdir $dirnew
	cp ${dir}/run.sbatch ${dirnew}
	#cp ${dir}/param.txt ${dirnew}
	RewriteParams ${dir}/param.txt ${dirnew}/param.txt
	if [ -e ${dir}/saclistR.txt ]; then
		cp ${dir}/saclistR.txt ${dirnew}
		cp ${dir}/saclistL.txt ${dirnew}
	fi
	cp ${dir}/station_*_*_*.txt ${dirnew}
	mv $dir ${dir}"~"
	echo $dirnew $dir >> $list_rerun
	sbatch ${dirnew}/run.sbatch
done
echo $list_rerun
