#!/bin/bash

if [ $# != 2 ]; then
	echo "Usage: "$0" [PosteriorD file] [Emax (search for minimum when Emax<=0)]"
	exit
fi

fin=$1
Emax=$2

if [ `echo $Emax | awk '{if($1<=0.){print 0}else{print 1}}'` == 1 ]; then
	grep "search\#" $fin | grep -v 'rej' | awk 'NF>=20' | sed s/'(accepted)'/''/ | sed s/'(best)'/''/ | awk -v Emax=$Emax -F"E = " '{if($2<Emax)print $0}'
else
	grep "search\#" $fin | grep -v 'rej' | awk 'NF>=20'| sed s/'(accepted)'/''/ | sed s/'(best)'/''/ | awk -F"E = " 'BEGIN{Emin=1e30;line=NaN}{if($2<Emin){Emin=$2; line=$0}}END{print line}'
fi
#tail -n2 PosteriorD.txt; echo ""
