#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -c 1# specify 6 threads per process
#SBATCH -t 00:30:00
#SBATCH -p workq
#SBATCH -A loni_proteus01s
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err
# below are job commands
date

#export WORK_DIR=/work/$USER
module purge
#module load proteus/1.7.5
module load proteus/1.8.1

mkdir -p $WORK/floating_body_H.$SLURM_JOBID
cd $WORK/floating_body_H.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/floating_body_H.py .
cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
cp $SLURM_SUBMIT_DIR/floating_body_H.sh .

srun -N1 -n48 -c1 parun --TwoPhaseFlow floating_body_H.py -l 5 -C "he=0.01 T=30. fr=1.0" -O petsc.options.superlu_dist

date

exit 0

