#!/bin/bash

ProduceSBATCH() {
local _fparam=$1
local _logdir=$2
local _Rtime=$3
local _outname=$4
local _addoption=$5
local _nthreads=12
#local _run_exe=/projects/yeti4009/code/Programs/ExploitEvent/SearchLocation/SearchSource_SA_MC
local _run_exe=/projects/yeti4009/eqkhyposolver/EQKSolver
#local _run_exe=/home/tianye/EQKHypoSolver/EQKSolver
echo "#!/bin/bash

#SBATCH -J EQKSearch_${_logdir}
#SBATCH -o ${_logdir}/log_EQKSearch_%j.out
#SBATCH -e ${_logdir}/log_EQKSearch_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=${_nthreads}
#SBATCH --time=${_Rtime}:00:00
#SBATCH --mem=1gb

ulimit -c unlimited
module load gcc/gcc-4.9.2

export OMP_NUM_THREADS=${_nthreads}; time ${_run_exe} ${_fparam} ${_addoption}
" > ${_outname}

#export OMP_NUM_THREADS=${_nthreads}; time ${_run_exe} ${_fparam} ${_addoption} 1> ${_logdir}/log_EQKSearch.out 2> ${_logdir}/log_EQKSearch.err
#nohup bash -c \"export OMP_NUM_THREADS=${_nthreads}; (time ${_run_exe} ${_fparam}) 1> ${_logdir}/log.out 2> ${_logdir}/log.err\" &
}


### main ###
if [ $# != 4 ]; then
	echo "Usage: "$0" [clon] [clat] [#sta] [#trial]"
	exit
fi

clon=$1
clat=$2
Nsta=$3
Ntrl=$4


for wtype in B R L; do # wave type ( B R L )
	for dtype in 1000; do # data type ( 500 1000 1000sparse 1000PIazi )
		for mtype in Ei 1D; do # model type ( Ei 3D 1D )
			for ((itrl=1; itrl<=$Ntrl; itrl++)); do
			### label
			_label=${wtype}_${dtype}_${mtype}_Sparse${itrl}
			### path to the vel maps
			if [ $mtype == "Ei" ]; then
				_mdir=VelMaps_Eikonal
				_msuf=txt
				fparam=param_base.txt
			elif [ $mtype == "3D" ]; then
				_mdir=VelMaps_Model
				_msuf=txt
				fparam=param_base.txt
			elif [ $mtype == "1D" ]; then
				_mdir=VelMaps_Eikonal
				_msuf=txt_avg
				fparam=param_base_1D.txt
			else
				echo "Unknown model type: "$mtype; exit
			fi
			fsta=`/projects/yeti4009/eqkhyposolver/Scripts/PickStaSparse.sh $clon $clat $Nsta`
			### running time
			Rtime=`echo $Nsta | awk '{nhr=0.5+$1*0.1; if(nhr>12){nhr=12} printf "%.0f\n", nhr}'`
			### check param file
			if [ ! -e $fparam ]; then
				echo "no "$fparam" found!"
				exit
			fi
			# directory
			_dir=results_SAMC_${_label}
			if [ -e ${_dir} ]; then
				mkdir -p old_results/${_dir}
				ls ${_dir}/* | xargs -I file mv file old_results/${_dir}
			fi
			mkdir -p ${_dir}
			# station list (sparse)
			mv $fsta $_dir
			fsta=${_dir}/${fsta}
			# produce fparam
			#_dis_measure_label=`echo $dis | sed s/'_freeF'/''/`
			more $fparam  | sed s/'dflag base'/'dflag '${wtype}/ | sed s/'VelMaps'/${_mdir}/g | sed s/'txt_base'/${_msuf}/g | sed s/'_dis500'/'_dis'${dtype}/g | sed s/'results_SAMC_default'/${_dir}/g | awk -v fsta=$fsta '{if($1=="fRm"||$1=="fLm"){print $0,fsta}else{print $0}}' > ${_dir}/param.txt
			ProduceSBATCH ${_dir}/param.txt ${_dir} $Rtime ${_dir}/run.sbatch -c
			# run/submit
			echo Starting ${_dir}...
			sbatch ${_dir}/run.sbatch
			#sh ${_dir}/run.sbatch
			sleep 1
			done # itrl
		done
	done
done

echo "Waiting for results..."

while [ 1 == 1 ]; do
	if [ `squeue -u yeti4009 | wc -l` -le 1 ]; then break; fi
	sleep 60
done
dir_sparse=results_Sparse
if [ -e ${dir_sparse} ]; then
	rm -rf old_results/${dir_sparse}
	mv ${dir_sparse} old_results
fi
mkdir ${dir_sparse}
ls -d results_SAMC_?_*_*_Sparse* | xargs -I dir mv dir ${dir_sparse}

echo "All done!"

