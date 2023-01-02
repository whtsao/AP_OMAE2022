#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -c 1# specify 6 threads per process
#SBATCH -t 72:00:00
#SBATCH -p workq
#SBATCH -A loni_proteus01r
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err
# below are job commands
date

export WORK_DIR=/work/$USER

module purge
module load proteus/1.7.5

mkdir -p $WORK_DIR/floating_body_H.$SLURM_JOBID
cd $WORK_DIR/floating_body_H.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/floating_body_H.py .
cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
cp $SLURM_SUBMIT_DIR/floating_body_H.sh .

srun -N1 -n48 -c1 parun --TwoPhaseFlow floating_body_H.py -l 5 -C "he=0.01 T=30. fr=0.9" -O petsc.options.superlu_dist

date

exit 0

