#!/bin/bash -l
# Fine time sampling for the time slider: one array task per snapshot (extract_movie.py)
#   sbatch run_movie.sbatch            # snapshots 27 .. 139
#   sbatch --array=50,51 run_movie.sbatch
#SBATCH -o /ptmp/uli/coma_cubes/movie/slurm_%A_%a.out
#SBATCH -e /ptmp/uli/coma_cubes/movie/slurm_%A_%a.out
#SBATCH -D ./
#SBATCH -J movie
#SBATCH --array=27-139
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=90000
#SBATCH --time=00:40:00
module purge
module load anaconda/3/2023.03
export OMP_NUM_THREADS=1
export SLURM_CPUS_PER_TASK=16
python extract_movie.py ${SLURM_ARRAY_TASK_ID}
