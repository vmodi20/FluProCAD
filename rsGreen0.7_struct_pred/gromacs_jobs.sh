#!/bin/bash
#SBATCH --job-name=FPC-mut
#SBATCH --account=project_2001312
#SBATCH --partition=large
#SBATCH --time=24:00:00
#SBATCH --nodes=4
#SBATCH --tasks-per-node=40
#SBATCH --mem-per-cpu=4000

# this script runs a 192 core (8 full nodes) gromacs job.
export OMP_NUM_THREADS=1

module load gcc/9.1.0
module load hpcx-mpi/2.4.0
module load gromacs/2020.3       #change if you want a different version

#gmx_mpi grompp -f em.mdp -p topol.top -c ions.gro -o em
#srun gmx_mpi mdrun -deffnm em

#gmx_mpi grompp -f nvt.mdp -p topol.top -c em.gro -r em.gro -o nvt
#srun gmx_mpi mdrun -deffnm nvt

gmx_mpi grompp -f npt.mdp -p topol.top -c nvt.gro -t nvt.cpt -r nvt.gro -o npt
srun gmx_mpi mdrun -deffnm npt

gmx_mpi grompp -f md.mdp -p topol.top -c npt.gro -t npt.cpt -o topol
srun gmx_mpi mdrun -s topol -dlb yes
