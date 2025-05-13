#!/bin/bash
#PBS -l select=2:ncpus=16
#PBS -N PREDICT
#PBS -q XXX
#PBS -V

# === initiate ===
cd /home/hpc/chiayicheng/scratch/users/Kai/PREDICT  ### Your "PREDICT" folder path
source /home/hpc/chiayicheng/package/miniforge3/etc/profile.d/conda.sh
conda activate PREDICT

# === Example for each module ===
# For the meaning of each argument, please see README.md "Arguments".

# --- FindKmers ---
# Basic usage (use default arguments)
python /home/hpc/chiayicheng/scratch/users/Kai/PREDICT/PREDICT.py FindKmers

# Custom arguments
# python /home/hpc/chiayicheng/scratch/users/Kai/PREDICT/PREDICT.py FindKmers --features gene --pThreshold 0.01 --alg RandomForest --up_stream 2000 --down_stream 800


# --- Kmer2Motif ---
# Basic usage (use default arguments)
# python /home/hpc/chiayicheng/scratch/users/Kai/PREDICT/PREDICT.py Kmer2Motif

# Custom arguments
# python /home/hpc/chiayicheng/scratch/users/Kai/PREDICT/PREDICT.py Kmer2Motif --TopKmers 0.05 --KeepTopMotifs 0.2 --ScoreCutOff 0.85


# --- RanKmers ---
# Basic usage (use default arguments)
# python /home/hpc/chiayicheng/scratch/users/Kai/PREDICT/PREDICT.py RanKmers

# Custom arguments
# python /home/hpc/chiayicheng/scratch/users/Kai/PREDICT/PREDICT.py RanKmers --features gene --alg RandomForest --up_stream 300 --down_stream 700


# --- TeamKmers ---
# Basic usage (use default arguments)
# python /home/hpc/chiayicheng/scratch/users/Kai/PREDICT/PREDICT.py TeamKmers

# Custom arguments
# python /home/hpc/chiayicheng/scratch/users/Kai/PREDICT/PREDICT.py TeamKmers --Motif False --TPCSThreshold 0.75 --TPTNCSThreshold 


# --- ViewKmers --- NOT recommended to run this module on HPC, better to run it locally in terminal.
# Basic usage (use default arguments)
# streamlit run /home/hpc/chiayicheng/scratch/users/Kai/PREDICT/PREDICT.py -- ViewKmers

# Custom arguments
# streamlit run /home/hpc/chiayicheng/scratch/users/Kai/PREDICT/PREDICT.py -- ViewKmers --features gene --up_stream 1000 --down_stream 500



