#!/bin/bash

# file lists to link
filelist1=(src_SDContainer/Point.h src_SDContainer/Map.h src_SDContainer/Array2D.h src_SDContainer/DisAzi.h src_SDContainer/Map.cpp src_SDContainer/StaList.h src_MISC/MyOMP.h src_MISC/StackTrace.h src_RadPattern/RadPattern.h src_RadPattern/RadPattern.cpp)

filelist2=(/projects/yeti4009/code/MyLibs/include/Point.h /projects/yeti4009/code/MyLibs/include/Map.h /projects/yeti4009/code/MyLibs/include/Array2D.h /projects/yeti4009/code/MyLibs/include/DisAzi.h /projects/yeti4009/code/MyLibs/sources/Map/Map.cpp /projects/yeti4009/code/MyLibs/include/StaList.h /projects/yeti4009/code/MyLibs/include/MyOMP.h /projects/yeti4009/code/MyLibs/include/StackTrace.h /projects/yeti4009/code/Programs/RadPattern/RadPattern_C/RadPattern.h /projects/yeti4009/code/Programs/RadPattern/RadPattern_C/RadPattern.cpp)

ifile=${#filelist1[@]}
for file1 in `ls ./src_RadPattern/*.f`; do
	file2=`echo $file1 | awk -F/ '{print "/projects/yeti4009/code/Programs/RadPattern/RadPattern_C/RadPattern_src/"$NF}'`
	filelist1[ifile]=$file1
	filelist2[ifile]=$file2
	let ifile++
done


# check differences
Nlinediff=0
for ((ifile=0;ifile<${#filelist1[@]};ifile++)); do
	file1=${filelist1[ifile]}; file2=${filelist2[ifile]}
	if [ ! -e $file1 ] || [ ! -e $file2 ]; then
		echo "file "$file1" or "$file2" not found!"
		exit
	fi
	Nlinediff=$((Nlinediff + `sdiff -s $file1 $file2 | wc -l`))
	echo -e "Checked: "$file1 - $file2"\n  : total Ndiff = "$Nlinediff
done

if [ $Nlinediff -gt 0 ]; then
	echo "Stopped: difference(s) detected!!"
	exit
fi

# link
mkdir -p trash
for ((ifile=0;ifile<${#filelist1[@]};ifile++)); do
	mv ${filelist1[ifile]} trash/
	#rm -f ${filelist1[ifile]}
	ln ${filelist2[ifile]} ${filelist1[ifile]}
done

echo "All linked! old files moved to ./trash!"
