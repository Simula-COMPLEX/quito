#!/bin/bash
#SBATCH -p dgx2q
#SBATCH --job-name=Quito
#SBATCH -n 10
#SBATCH --time 00-48:00:00    # time (D-HH:MM:SS)
module purge
module load slurm/20.02.7
source /home/xinyi/anaconda3/bin/activate quito_gpu
export PYTHONPATH=root_path:$PYTHONPATH
python /home/xinyi/quito/ex3/quito_test.py
