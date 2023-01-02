#!/bin/bash
#SBATCH -N 4
#SBATCH -n 192
#SBATCH -c 1 # specify 6 threads per process
#SBATCH -t 00:30:00
#SBATCH -p workq
#SBATCH -A loni_proteus01s
#SBATCH -o o.out # optional, name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e e.err
# below are job commands
date

#export WORKDIR=/work/$USER

module purge
module load proteus/1.7.5

mkdir -p $WORKDIR/P_uncon_floatbody.$SLURM_JOBID
cd $WORKDIR/P_uncon_floatbody.$SLURM_JOBID 
cp $SLURM_SUBMIT_DIR/*.py .
cp $SLURM_SUBMIT_DIR/petsc.options.superlu_dist .
cp $SLURM_SUBMIT_DIR/petsc.options.asm .
cp $SLURM_SUBMIT_DIR/*.sh .

srun -N4 -n192 -c1 parun --TwoPhaseFlow floating_body.py -l 5 -C "final_time=20.0 he=0.005 fr=1.0" -O petsc.options.superlu_dist
#srun parun --TwoPhaseFlow floating_body.py -F -l 5 -C "he=0.002 fr=1.0" -O petsc.options.asm

date

exit 0

