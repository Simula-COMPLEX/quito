#!/bin/bash
#SBATCH -p defq
#SBATCH --job-name=Quito
#SBATCH -n 10
#SBATCH --time 00-48:00:00    # time (D-HH:MM:SS)
module purge
module load slurm/20.02.7
source /home/xinyi/anaconda3/bin/activate quantum_env
export PYTHONPATH=root_path:$PYTHONPATH
python ../Quito_CoverageRunning/quito.py