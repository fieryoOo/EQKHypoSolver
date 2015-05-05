#!/bin/bash

ProduceSBATCH() {
local _fparam=$1
local _logdir=$2
local _Rtime=$3
local _outname=$4
local _addoption=$5
local _nthreads=12
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


### main ###
for type in B R L; do # R L B C
   for dis in 1000; do # 1000 500 1000_sparse 1000_PIazi 1000_freeF
      for mtype in 3d avg; do #3d avg
			### label
			if [ $mtype == "3d" ]; then
				_label=${type}${dis}
				_ftag=txt
			else
				_label=${type}${dis}_avg
				_ftag=txt_avg
			fi
			### time
			if [ ${_label} == "L500" ] || [ ${_label} == "R500" ]; then
				Rtime=5
			else
				Rtime=15
			fi
			### param file
			fparam=param_base.txt
			if [ ! -e $fparam ]; then
				echo "no param_base.txt found!"
				exit
			fi
			for _VelMapType in Eikonal; do # Model Eikonal
				if [ $_VelMapType == "Eikonal" ]; then
					_label=${_label}_Eikonal
				fi
				# directory
				_dir=results_SAMC_${_label}
				if [ -e ${_dir} ]; then
					mkdir -p old_results/${_dir}
					ls ${_dir}/* | xargs -I file mv file old_results/${_dir}
				fi
				mkdir -p ${_dir}
				# fparam
				_dis_measure_label=`echo $dis | sed s/'_freeF'/''/`
				more $fparam  | sed s/'dflag base'/'dflag '${type}/ | sed s/'VelMaps'/'VelMaps_'${_VelMapType}/g | sed s/'txt_base'/${_ftag}/g | sed s/'_dis500'/'_dis'${_dis_measure_label}/g | sed s/'results_SAMC_default'/${_dir}/g > ${_dir}/param.txt
				ProduceSBATCH ${_dir}/param.txt ${_dir} $Rtime ${_dir}/run.sbatch -c
				echo Starting ${_dir}...
				#sbatch ${_dir}/run.sbatch
				sh ${_dir}/run.sbatch
				sleep 1
			done
      done
   done
done
