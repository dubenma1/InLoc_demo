#!/bin/bash
#SBATCH --job-name=InLoc_demo
#SBATCH --nodes=1
#SBATCH --partition gpu
#SBATCH --output=InLoc_demo.log
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --time=1-00:00:00
module load MATLAB/2018a
nvidia-smi
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/lucivpav/gflags-2.2.2/build/lib
cat startup.m inloc_demo.m | matlab -nodesktop
