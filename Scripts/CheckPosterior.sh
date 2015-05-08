#!/bin/bash

if [ $# != 1 ]; then
	echo "Usage: "$0" [Emax]"
	exit
fi

grep -v 'rej' PosteriorD.txt  | sed s/'(accepted)'/''/ | sed s/'(best)'/''/ | awk -v Emax=$1 -F"E = " '{if($2<Emax)print $0}'
#tail -n2 PosteriorD.txt; echo ""
