#!/bin/bash

# file lists to link
filelist1=(./Point.h ./Map.h ./Array2D.h ./DisAzi.h ./MyOMP.h ./StackTrace.h ./Map.cpp ./RadPattern.h ./RadPattern.cpp ./StaList.h)
filelist2=(/projects/yeti4009/code/MyLibs/include/Point.h /projects/yeti4009/code/MyLibs/include/Map.h /projects/yeti4009/code/MyLibs/include/Array2D.h /projects/yeti4009/code/MyLibs/include/DisAzi.h /projects/yeti4009/code/MyLibs/include/MyOMP.h /projects/yeti4009/code/MyLibs/include/StackTrace.h /projects/yeti4009/code/MyLibs/sources/Map/Map.cpp /projects/yeti4009/code/Programs/RadPattern/RadPattern_C/RadPattern.h /projects/yeti4009/code/Programs/RadPattern/RadPattern_C/RadPattern.cpp /projects/yeti4009/code/MyLibs/include/StaList.h)
ifile=${#filelist1[@]}
for file1 in `ls ./RadPattern_src/*.f`; do
	file2=/projects/yeti4009/code/Programs/RadPattern/RadPattern_C/${file1}
	filelist1[ifile]=$file1
	filelist2[ifile]=$file2
	let ifile++
done


# check differences
Nlinediff=0
for ((ifile=0;ifile<${#filelist1[@]};ifile++)); do
	Nlinediff=$((Nlinediff + `sdiff -s ${filelist1[ifile]} ${filelist2[ifile]} | wc -l`))
	echo "Checked: "${filelist1[ifile]} ${filelist2[ifile]}": total Ndiff = "$Nlinediff
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
