#!/bin/bash -l
#SBATCH --job-name=bench_intranode
#SBATCH --output=Fritz_ICX_MultiGrid
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --cpu-freq=2400000-2400000:performance

unset SLURM_EXPORT_ENV

module load likwid intel intelmpi
export I_MPI_PIN=1
export I_MPI_DEBUG=0

FILENAME="Fritz_ICX_MultiGrid.csv"

cd ~/HighPerformanceComputingSolver
make distclean
make

rm $FILENAME
touch $FILENAME

function _iterate() {
    # Bash Array Definition
    declare -a process_counts=(1 2 4 8 15 30 40 50 61 72)
    for npn in "${process_counts[@]}"; do
        np_1=$(($npn - 1))
        export I_MPI_PIN_PROCESSOR_LIST=0-$np_1

        result="$(srun -n $npn ./exe-ICX dcavity.par)"
        result="$(echo $result | sed 's/MPI startup(): Warning: I_MPI_PMI_LIBRARY will be ignored since the hydra process manager was found //g')"
        
        echo $npn $result >>$FILENAME
    done
}

_iterate
