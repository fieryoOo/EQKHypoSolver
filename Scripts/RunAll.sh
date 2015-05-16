#!/bin/bash

ProduceSBATCH() {
local _fparam=$1
local _logdir=$2
local _Rtime=$3
local _outname=$4
local _addoption=$5
local _nthreads=8
#local _run_exe=/projects/yeti4009/code/Programs/ExploitEvent/SearchLocation/SearchSource_SA_MC
#local _run_exe=/projects/yeti4009/eqkhyposolver/EQKSolver
local _run_exe=/home/tianye/EQKHypoSolver/EQKSolver
echo "#!/bin/bash

#SBATCH -J EQKSearch_${_logdir}
#SBATCH -o ${_logdir}/log_EQKSearch_%j.out
#SBATCH -e ${_logdir}/log_EQKSearch_%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=${_nthreads}
#SBATCH --time=${_Rtime}:00:00
#SBATCH --mem=1gb

ulimit -c unlimited
module load gcc/4.8.0

export OMP_NUM_THREADS=${_nthreads}; time ${_run_exe} ${_fparam} ${_addoption} 1> ${_logdir}/log_EQKSearch.out 2> ${_logdir}/log_EQKSearch.err
" > ${_outname}

#nohup bash -c \"export OMP_NUM_THREADS=${_nthreads}; (time ${_run_exe} ${_fparam}) 1> ${_logdir}/log.out 2> ${_logdir}/log.err\" &
}


### main 
for wtype in B R; do # wave type ( B R L )
	for dtype in 1000 500; do # data type ( 500 1000 1000sparse 1000PIazi )
		for mtype in Ei 1D; do # model type ( Ei 3D 1D )
			### label
			_label=${wtype}_${dtype}_${mtype}
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
			### running time
			if [ ${_label} == "L500" ] || [ ${_label} == "R500" ]; then
				Rtime=5
			else
				Rtime=15
			fi
			### check param file
			if [ ! -e $fparam ]; then
				echo "no param_base.txt found!"
				exit
			fi
			# directory
			_dir=results_SAMC_${_label}
			if [ -e ${_dir} ]; then
				mkdir -p old_results/${_dir}
				ls ${_dir}/* | xargs -I file mv file old_results/${_dir}
			fi
			mkdir -p ${_dir}
			# produce fparam
			#_dis_measure_label=`echo $dis | sed s/'_freeF'/''/`
			more $fparam  | sed s/'dflag base'/'dflag '${wtype}/ | sed s/'VelMaps'/${_mdir}/g | sed s/'txt_base'/${_msuf}/g | sed s/'_dis500'/'_dis'${dtype}/g | sed s/'results_SAMC_default'/${_dir}/g > ${_dir}/param.txt
			ProduceSBATCH ${_dir}/param.txt ${_dir} $Rtime ${_dir}/run.sbatch -c
			# run/submit
			echo Starting ${_dir}...
			#sbatch ${_dir}/run.sbatch
			sh ${_dir}/run.sbatch
			sleep 1
		done
	done
done
