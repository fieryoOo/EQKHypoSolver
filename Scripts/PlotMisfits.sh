#!/bin/bash

PlotText() {
   local _REG=$1
   local _text=$2
   local _XS=$3
   local _YS=$4

   #local _lb=('a' 'b' 'c' 'd' 'e' 'f')
   #local _title=`echo ${_text} | awk -v lb=${_lb[ifile]} '{print "("lb")  "$0}'`
   local _title=$_text
   echo ${_title}
   local _llon=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v xs=$_XS '{print $1+xs*($2-$1)}'`
   local _ulat=`echo $_REG | sed s/'\-R'/''/ | awk -F/ -v ys=$_YS '{print $4+ys*($4-$3)}'`
   #echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "15. 0. 20 LT", $0}' | pstext -R -J -Wlightgray,O3 -O -K -N >> $psout
   echo $_title | awk -v llon=$_llon -v ulat=$_ulat '{print llon, ulat, "12. 0. 20 LT", $0}' | pstext -R -J -O -K -N >> $psout
   let ifile++
}

GetSNR() {
	local _sta=$1
	if [ $useW == true ]; then
		GetSNR1 $_sta
	else
		GetSNR2 $_sta
	fi
}

GetSNR1() {
	local _sta=$1
	local _fsac=`ls ${wdir}/${_sta}.?HZ_real.SAC 2> /dev/null`
	succ=true
	if [ `echo $_fsac | awk '{print NF}'` == 0 ] || [ ! -e $_fsac ]; then
		echo "sac file "$_fsac" not found for SNR"
		succ=false; return
	fi
	SNR=`~/usr/bin/saclst user1 f $_fsac | awk '{print $2}'`
	if [ $SNR == -12345 ]; then
		echo "failed to read SNR from sacheader.user1"
		succ=false; return
	fi
}

GetSNR2() {
	succ=true
   local _sta=$1
   local _f1=/work1/tianye/EQKLocation/SAC/Disps/${event}_disp_LHZ_10sec.txt
   local _f2=/work1/tianye/EQKLocation/SAC/Disps/${event}_disp_LHZ_16sec.txt
	if [ ! -e $_f1 ] || [ ! -e $_f2 ]; then
		_f1=/lustre/janus_scratch/yeti4009/EQKLocation/SAC/Disps/${event}_disp_LHZ_10sec.txt
		_f2=/lustre/janus_scratch/yeti4009/EQKLocation/SAC/Disps/${event}_disp_LHZ_16sec.txt
		if [ ! -e $_f1 ] || [ ! -e $_f2 ]; then succ=false; return; fi
	fi

	#if [ $per == 10 ]; then
	#	SNR=`awk -v sta=$_sta '$4==sta{print $11}' $_f1`
	#else
	#	SNR=`awk -v sta=$_sta '$4==sta{print $11}' $_f2`
	#fi
   #SNR=`echo $SNR | awk '{ if(NF<1){ snr=2 }else{ snr=$1; if(snr<2){snr=2} } print snr} '`

   SNR1=`awk -v sta=$_sta '$4==sta{print $11}' $_f1`
   SNR2=`awk -v sta=$_sta '$4==sta{print $11}' $_f2`
   SNR=`echo $SNR1 $SNR2 | awk '{ if(NF<2){ snr=2 }else{ snr=$1; if(snr>$2){snr=$2} if(snr<2){snr=2} } print snr} '`
}

SNR2Size() {
	local _SNR=$1
	local _size=2
	if [ $useW == true ]; then
		_size=`echo $_SNR | awk '{size=0.008*$1; if(size>0.4){size=0.4;} print size}'`
	else
		_size=`echo $_SNR | awk '{size=0.025*$1; if(size>0.4){size=0.4;} print size}'`
	fi
	echo $_size
}

GetCC() {
	local _sta=$1
	local _fsac1=`ls ${wdir}/${_sta}.?HZ_real.SAC 2> /dev/null`
	local _fsac2=`ls ${wdir}/${_sta}.?HZ_syn.SAC 2> /dev/null`
	succ=true
	if [ ! -e $_fsac1 ] || [ ! -e $_fsac2 ]; then
		echo "sac file "$_fsac1" or "$_fsac2" not found for CC"
		succ=false; return
	fi
	CC=`$exeCC $_fsac1 $_fsac2`
}

PlotHisto() {
	local _fin=$1
	local _REG=$2
	local _colorflag=$3
	local _ptext=$4
	local _YS=$5
	pshistogram $_fin $_REG -J ${_colorflag} -Z1 -O -K >> $psout
	if [ $_ptext != true ]; then return; fi
	local _avg=`awk 'BEGIN{a=0}{a+=$1}END{print a/NR}' $_fin`
	local _std=`awk -v avg=$_avg 'BEGIN{std=0}{std+=($1-avg)**2}END{print (std/(NR-1))**0.5}' $_fin`
	local _label=`echo $_avg $_std | awk '{printf "%.2f +/- %.2f\n",$1,$2}'`
	PlotText $_REG "$_label" 0.1 ${_YS}
}

### compute/plot surface
PlotMap() {
	# input params
	local _fin=$1;	
	local _icol=$2
	local _label=$3
	# trash bin
	local ftobedeleted=( '' ); local idel=0
	# grab data
	local _fdata=${_fin}.tmp
	awk -v ic=${_icol} '{print $1,$2,$ic,$6}' ${_fin} > ${_fdata}
	ftobedeleted[idel]=$_fdata; let idel++
	# surface params
	local res=0.1; local ts=0.2; local bdis=30
	local fbmtmp1=${_fdata}.bm1
	ftobedeleted[idel]=$fbmtmp1; let idel++
	local fbmtmp2=${_fdata}.bm2
	ftobedeleted[idel]=$fbmtmp2; let idel++
	# redistribute data point locations
	blockmean ${_fdata} $REG -I$bdis'km' > $fbmtmp1
	blockmean $fbmtmp1 $REG -F -I$bdis'km' > $fbmtmp2
	# compute surface
	local fgrd=${_fdata}.grd
	ftobedeleted[idel]=$fgrd; let idel++
	surface $fbmtmp2 -T$ts -G$fgrd -I$res $REG
	# resample and plot
	local sms=`echo $res | awk '{print $1/5.0}'`
	local fgrds=${_fdata}.grds
	ftobedeleted[idel]=$fgrds; let idel++
	grdsample $fgrd -Q -G$fgrds $REG -I$sms
	grdimage -R -J $fgrds -C$fcpt -O -K >> $psout
	psbasemap -R -J -Bwesn -O -K >> $psout
	# boundaries
	pscoast -R -J -A100 -N1/3/0/0/0 -N2/3/0/0/0 -O -K -W3 >> $psout
	psxy ${Hdir}/wus_province_II.dat -R -J -W5/255/0/0 -M"99999 99999"  -O -K >> $psout
	psxy ${Hdir}/platebound.gmt -R -J -W5/255/0/0 -M"99999 99999"  -O -K >> $psout
	# station location and color bar
	#psxy ${_fdata} -R -J -Sc -W2,black -C$fcpt -O -K >> $psout
	sort -g -k3 ${_fdata} | psxy -R -J -Sc -W2,black -C$fcpt -O -K >> $psout
	#psscale  -C$fcpt -B4.5:"Misfit (sec)": -P -D3.65/-1./7.3/0.3h -O -K >> $psout
	psscale  -C$fcpt -B:"${_label}": -P -D3.6/-1./7.2/0.3h -O -K >> $psout
	# centre (earthquake) location
	clon=-114.9; clat=41.14
	echo $clon $clat | psxy -R -J -Sa0.6 -W2,blue -Gred -O -K >> $psout
	# SNR-size legend
	pslegend -R -D-123.3/32.2/2.4/2.2/BL -J -F -G220 -O -K <<- EOF >> $psout
	G 0.05i
	S 0.1i c ${sizes[0]} white 0p 0.3i SNR=${SNRs[0]}
	G 0.1i
	S 0.1i c ${sizes[1]} white 0p 0.3i SNR=${SNRs[1]}
	G 0.1i
	S 0.1i c ${sizes[2]} white 0p 0.3i SNR=${SNRs[2]}
	G 0.05i
	EOF
	# clean-up
	echo ${ftobedeleted[@]} | xargs rm -f
}


### main ###
if [ $# != 1 ] && [ $# != 2 ] && [ $# != 3 ]; then
	echo "Usage: "$0" [result_file] [waveform dir (optional)] [plot type (default=0=group, 1=CorrC)]"
	exit
fi

fin=$1
if [ ! -e $fin ]; then
	echo "bad file: "$fin
	exit
fi

# event
event=`pwd | awk -F/ '{print $(NF-1)}'`

# period
#per=`echo $fin | sed s/'_'/' '/g | sed s/'\.'/' '/g | awk '{for(i=1;i<10;i++)print $i}' | grep sec | sed s/'sec'/''/`

# file label
#label=`readlink -f $fin | awk -F/ '{print $(NF-1)"_"$(NF)}' | cut -d. -f1 | awk -F_ '{print "data-used="$3" Map="$5" data-plot="$6" "$10}'`
label=`echo $fin | cut -d. -f1 | awk -F_ '{print $1,$(NF)}'`
echo "label = "$label

### set up gmt
gmtset HEADER_FONT_SIZE 12
gmtset LABEL_FONT_SIZE 12
gmtset ANNOT_FONT_SIZE 10
gmtset HEADER_OFFSET 0.

### find script/header directory
Sdir=/home/tianye/EQKHypoSolver/Scripts
if [ ! -e $Sdir ]; then
	Sdir=/projects/yeti4009/eqkhyposolver/Scripts
	if [ ! -e $Sdir ]; then echo "Sdir not found"; exit; fi
fi

Hdir=/home/tianye/code/Programs/head
if [ ! -e $Hdir ]; then
	Hdir=/projects/yeti4009/code/Programs/head
	if [ ! -e $Hdir ]; then echo "Hdir not found"; exit; fi
fi

exeCC=`which SAC_Correlation`
if [ ! -e $exeCC ]; then
	echo "Cannot find SAC_Correlation"
	exit
fi

fsta=${Sdir}/station.lst

### check existence of dir ./waveforms
if [ $# -ge 2 ]; then
	useW=true
	wdir=$2
elif [ -e ./waveforms ]; then
	useW=true
	wdir=./waveforms
	echo "results recognized as Waveform fitting"
elif [ -e ./waveforms_init ]; then
	useW=true
	wdir=./waveforms_init
	echo "results recognized as Waveform fitting"
else
	useW=false
	echo "results recognized as FTAN measurements"
fi

### compute
if [ $# == 3 ] && [ $3 == 1 ]; then
	cc=true
else
	cc=false
fi
Ntmp=`grep -n '#' $fin | awk -F: '{print $1}' | tail -n1`
#awk -v N=$Ntmp '{if(NR>N){print $1,$2,$8-$9-$10}}' $fin | awk 'NF==3' | psxy -R -J -Sc0.5 -W1,black -C$fcpt -O -K >> $psout
fdata=${fin}.tmp; rm -f $fdata
awk -v N=$Ntmp '{if(NR>N&&$9>-999&&$10>-999){print $1,$2,$5-$6-$7,$8-$9-$10,($11-$12)/$12,$3,$4}}' $fin | awk 'NF==7' | while read lon lat grTmis phTmis ampPerc dis azi; do
#awk -v N=$Ntmp '{if(NR>N&&$9>-999&&$10>-999){print $1,$2,$5-$6,$7}}' $fin | awk 'NF==4' | while read lon lat Tmis Tsrc; do
	sta=`awk -v lon=$lon -v lat=$lat 'BEGIN{if(lon<0){lon+=360.}}{lonsta=$2; if(lonsta<0.){lonsta+=360.} latsta=$3; if( (lon-lonsta)**2+(lat-latsta)**2 < 0.0001 ){print $1} }' $fsta`
	GetSNR $sta
	if [ $succ == false ]; then exit; fi
	size=`SNR2Size $SNR`
	if [ $cc == true ]; then
		GetCC $sta
		if [ $succ == false ]; then exit; fi
		grTmis=$CC
	fi
	echo $lon $lat $grTmis $phTmis $ampPerc $size $SNR $sta $dis $azi >> $fdata
done

# station sizes (SNR - size)
SNRs=( 0 0 0 )
if [ $useW == true ]; then
	SNRs=( 15 30 50 )
else
	SNRs=( 3 8 15 )
fi
sizes=( '' '' '' )
for ((i=0;i<3;i++)); do
   #sizes[i]=`echo ${SNRs[i]} | awk '{size=0.01*$1; if(size>0.5){size=0.5;} print size}'`
	sizes[i]=`SNR2Size ${SNRs[i]}`
done
echo "snrs: "${SNRs[@]}
echo "size: "${sizes[@]} 

### plot
psout=${fin}.ps
REG=-R-123/-108/32.5/46.5
SCA=-Jn-115.5/0.5

### group T
if [ $cc == true ]; then
	fcpt=${Sdir}/CC.cpt
	psbasemap $REG $SCA -Ba5f1/a5f1:."Correlation Coef ($label)":WeSn -X4.5 -Y6 -K > $psout
	cptlabel="Correlation"
else
	psbasemap $REG $SCA -Ba5f1/a5f1:."Group Time Misfit ($label)":WeSn -X4.5 -Y6 -K > $psout
	cptlabel="Misfit (sec)"
	if [ $useW == true ]; then
		fcpt=${Sdir}/GMisfit.cpt
	else
		fcpt=${Sdir}/GMisfit.cpt
	fi
fi
PlotMap $fdata 3 ${cptlabel}


### phase T
if [ $useW == true ]; then
	fcpt=${Sdir}/PMisfitW.cpt
	psbasemap -R -J -Ba5f1/a5f1:."Phase Misfit ($label)":WeSn -X8 -K -O >> $psout
	unit=radian
else
	fcpt=${Sdir}/PMisfit.cpt
	psbasemap -R -J -Ba5f1/a5f1:."Phase Time Misfit ($label)":WeSn -X8 -K -O >> $psout
	unit=sec
fi
PlotMap $fdata 4 "Misfit ("${unit}")"
#psscale  -C$fcpt -B1.0:"Misfit (${unit})": -P -D3.65/-1./7.3/0.3h -O -K >> $psout


### amplitude misfits in percentage
if [ $useW == true ]; then
	fcpt=${Sdir}/AMisfitW.cpt
else
	fcpt=${Sdir}/AMisfit.cpt
fi
psbasemap -R -J -Ba5f1/a5f1:."Amplitude Misfit Percentage ($label)":WeSn -X8 -K -O >> $psout
PlotMap $fdata 5 "Misfit (%)"
#psscale  -C$fcpt -B0.5:"Misfit (%)": -P -D3.65/-1./7.3/0.3h -O -K >> $psout


# plot histogram
#REG=-R-8/8/0/40
#ftmp_histo=ftmp_histo.txt
#psbasemap $REG -JX6/8 -Ba2.f0.5/a10f2:."Group/Phase Time Shift (sec)":WeSn -X9 -O -K >> $psout
#awk '{if($7>15)print $4}' $fdata > $ftmp_histo
#PlotHisto $ftmp_histo $REG "-W0.4 -Gsteelblue -L3,steelblue3" true -0.05
#awk '{if($7>15)print $3}' $fdata > $ftmp_histo
#PlotHisto $ftmp_histo $REG "-W0.5 -L3,lightred" true -0.15

##awk '{if($6>10)print $3-$4}' $fdata | pshistogram -R-8/8/0/40 -JX6/8 -W0.4 -Gsteelblue -L3,steelblue3 -Z1 -X9 -O -K >> $psout
#awk '{if($6>10)print $3}' $fdata | pshistogram -R -J -W0.5 -Ba1.f0.2/a10f2:."Phase Time Shift (sec)":WeSn -L3,lightred -Z1 -O -K >> $psout
#rm -f $fdata
echo $fdata

### finalize
pwd | psxy -R -J -O >> $psout
echo $psout

