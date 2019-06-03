#!/bin/tcsh
#SBATCH -J plot_aero
#SBATCH -n 1
#SBATCH --mem=50G
#SBATCH -t 24:00:00
#SBATCH -A P93300665
#SBATCH -p dav
#SBATCH -e err.%J
#SBATCH -o out.%J

### Initialize the Slurm environment
setenv  TMPDIR /glade/scratch/duvivier/temp
mkdir -p $TMPDIR

module load ncl/6.4.0

./plot_cesm2_atm_aero_spatial.csh

