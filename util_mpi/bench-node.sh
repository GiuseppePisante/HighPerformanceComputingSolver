#!/bin/bash -l
#SBATCH --job-name=bench_intranode
#SBATCH --output=Fritz_ICX_DMVM_intranode
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

FILENAME="result_bench_intranode.csv"

cd ~/Pampi/ex4/skeleton
make distclean
make

rm $FILENAME
touch $FILENAME
echo "Ranks,IT,N,Time,Omega,Res" >>$FILENAME

_iterate() {
    # Bash Array Definition
    declare -a process_counts=(1 2 4 8 15 30 40 50 61 72)
    for npn in "${process_counts[@]}"; do
        np_1=$(($npn - 1))
        export I_MPI_PIN_PROCESSOR_LIST=0-$np_1

        result="$(mpirun -n $npn ./exe-ICC poisson.par)"
        result="$(echo $result | sed 's/MPI startup(): Warning: I_MPI_PMI_LIBRARY will be ignored since the hydra process manager was found //g')"
        
        echo $npn $result >>$FILENAME
    done
}


# Array mit omega Werten
declare -a omega_values=(1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1)

# Für jeden omega Wert
for omega in "${omega_values[@]}"; do
    # Ändere omega in der poisson.par Datei
    sed -i "s/omg.*[0-9.]\+/omg      $omega/" poisson.par
    
    # Führe Iteration aus
    _iterate
done

sed -i 's/ /,/g' $FILENAME