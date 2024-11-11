#!/bin/bash
  
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -l ncpus=48
#PBS -l mem=190GB
#PBS -l jobfs=400GB
#PBS -l software=matlab_<institution>
#PBS -l wd
  
# Load modules, always specify version number.
module load matlab/R2019b
module load matlab_licence/<institution>
  
# Must include `#PBS -l storage=scratch/ab12+gdata/yz98` if the job
# needs access to `/scratch/ab12/` and `/g/data/yz98/`
  
# Run matlab application
matlab -nodisplay -nosplash -r "outputDir='$PBS_JOBFS',numberOfWorkers=$PBS_NCPUS, mfile, quit" > $PBS_JOBID.log