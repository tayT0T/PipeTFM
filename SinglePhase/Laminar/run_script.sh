#!/bin/bash
#SBATCH --job-name="Restart_bitch"
#SBATCH --output="slurm.%j.out"
#SBATCH --error="slurm.%j.err"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=0
#SBATCH --account="csd999"
#SBATCH --export=ALL
#SBATCH -t 4:30:00

module purge
module load cpu
#Load module file(s) into the shell environment
module load gcc/10.2.0
#module load openmpi/4.1.1
module load slurm

export OMPI_MCA_coll_hcoll_enable=0

# grid number, axial length, forcing, liquid height, U_GS, Froude, Reynold
mpirun -np 32 ./main_exe 512 4.5 0. 1.01 1.10 10.0 40.0
./ConvertingPPM.sh
./Combine.sh
