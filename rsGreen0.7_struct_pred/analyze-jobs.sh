#!/bin/bash
#SBATCH --job-name=FPC-analyze
#SBATCH --account=project_2001312
#SBATCH --partition=large
#SBATCH --time=24:10:00
#SBATCH --nodes=2
#SBATCH --tasks-per-node=40
#SBATCH --mem-per-cpu=4000

# this script runs a 192 core (8 full nodes) gromacs job.
export OMP_NUM_THREADS=1

module load gcc/9.1.0
module load hpcx-mpi/2.4.0
module load gromacs/2020.3       #change if you want a different version

#outtag=$1

./analyze-traj.sh $1 << EOF
5
3
TYG
0.30
EOF
