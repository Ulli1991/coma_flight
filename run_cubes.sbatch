#!/bin/bash -l
# Cube extraction for the Coma Virtual Observatory, one snapshot per job:
#   for s in 139 27 56 77 109 121; do sbatch --export=SNAP=$s -J cubes_$s run_cubes.sbatch; done
#SBATCH -o /ptmp/uli/coma_cubes/slurm_%x_%j.out
#SBATCH -e /ptmp/uli/coma_cubes/slurm_%x_%j.out
#SBATCH -D ./
#SBATCH -J cubes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=72
#SBATCH --mem=240000
#SBATCH --time=03:00:00
module purge
module load anaconda/3/2023.03
export OMP_NUM_THREADS=1
export SLURM_CPUS_PER_TASK=16
python extract_cubes.py ${SNAP:?set SNAP}
