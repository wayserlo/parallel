#!/bin/bash

#PBS -l walltime=00:10:00,nodes=5:ppn=3
#PBS -N task3
#PBS -q batch

cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=$PBS_NUM_PPN
./monte_carlo.out
