#!/bin/bash

if [ $# != 1 ]; then
	echo "Usage: "$0" [result_file]"
	exit
fi

fin=$1
if [ ! -e $fin ]; then
	echo "bad file: "$fin
	exit
fi

### plot
psout=${fin}.ps
REG=-R-123/-109/34/48
psbasemap $REG -Jn-116/1 -Ba5f1/a5f1:."Phase Misfits":WeSn -X5 -Y8 -K -P > $psout
pscoast -R -J -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 >> $psout
psxy /home/tianye/code/Programs/head/wus_province_II.dat -R -J -W5/255/0/0 -M"99999 99999"  -O -K >> $psout
psxy /home/tianye/code/Programs/head/platebound.gmt -R -J -W5/255/0/0 -M"99999 99999"  -O -K >> $psout

fcpt=Tshift.cpt
awk '{print $1,$2,$8-$9-$10}' $fin | awk 'NF==3' | psxy -R -J -Sc0.5 -W1,black -C$fcpt -O -K >> $psout

psscale  -C$fcpt -B1:"Phase Shift (sec)": -P -D6./-1.5/12./0.5h -O -K >> $psout

### finalize
pwd | psxy -R -J -O >> $psout
echo $psout

