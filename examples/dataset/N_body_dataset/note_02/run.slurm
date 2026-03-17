#!/bin/bash
#SBATCH --job-name=note_02_long
#SBATCH --partition=GPUA
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4         
#SBATCH --nodelist=node70           # ✅ 指定一个空闲节点
#SBATCH --time=1200:00:00
#SBATCH --output=genga-test-%j.out

module load compiler/cuda/12.5

./gengamod -R -1
