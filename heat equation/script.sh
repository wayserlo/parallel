#!/bin/bash
#PBS -l walltime=00:010:00,nodes=7:ppn=4
#PBS -N astanova
#PBS -q batch
cd $PBS_O_WORKDIR
for N in 2000 10000 50000
do
for((i = 1; i < 13; i++))
do
mpirun --hostfile $PBS_NODEFILE -np $i ./teplo.out $N
echo ","
done
echo "-------"
done
