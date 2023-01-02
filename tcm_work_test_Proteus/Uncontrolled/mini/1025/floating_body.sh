#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 72:00:00
#SBATCH -p workq
#SBATCH -A loni_proteus01r
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err
# below are job commands
date

export WORKDIR=/work/$USER

module purge
module load proteus/1.7.5

mkdir -p $WORKDIR/P_floating_body.$SLURM_JOBID
cd $WORKDIR/P_floating_body.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/floating_body.py .
cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
cp $SLURM_SUBMIT_DIR/floating_body.sh .

srun -N1 -n48 -c1 parun --TwoPhaseFlow floating_body.py -l 5 -C "final_time=20.0 he=0.02 fr=1.025" -O petsc.options.superlu_dist

date

exit 0

