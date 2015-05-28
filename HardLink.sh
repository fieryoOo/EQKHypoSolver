#!/bin/bash

# check code directory
dirs=(${dir_code} /home/tianye/code)
dir_code='dir_NaN'
for dir in ${dirs[@]}; do
	if [ -e $dir ]; then
		dir_code=$dir
		break;
	fi
done
if [ $dir_code == "dir_NaN" ]; then
	echo "cannot locate the code directory!"
	exit
fi

# file lists to link
filelist1=(src_SDContainer/Point.h src_SDContainer/Map.h src_SDContainer/Array2D.h src_SDContainer/DisAzi.h src_SDContainer/Map.cpp src_SDContainer/StaList.h src_MISC/MyOMP.h src_MISC/StackTrace.h src_RadPattern/RadPattern.h src_RadPattern/RadPattern.cpp)

filelist2=(${dir_code}/MyLibs/include/Point.h ${dir_code}/MyLibs/include/Map.h ${dir_code}/MyLibs/include/Array2D.h ${dir_code}/MyLibs/include/DisAzi.h ${dir_code}/MyLibs/sources/Map/Map.cpp ${dir_code}/MyLibs/include/StaList.h ${dir_code}/MyLibs/include/MyOMP.h ${dir_code}/MyLibs/include/StackTrace.h ${dir_code}/Programs/RadPattern/RadPattern_C/RadPattern.h ${dir_code}/Programs/RadPattern/RadPattern_C/RadPattern.cpp)

ifile=${#filelist1[@]}
for file1 in `ls ./src_RadPattern/*.f`; do
	file2=`echo $file1 | awk -v dir=${dir_code} -F/ '{print dir"/Programs/RadPattern/RadPattern_C/RadPattern_src/"$NF}'`
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
