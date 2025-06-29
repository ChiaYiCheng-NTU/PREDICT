#!/bin/bash
#PBS -l select=2:ncpus=16
#PBS -N PREDICT
#PBS -q XXX
#PBS -V

# === initiate ===
echo "conda base path:$CONDA_BASE"
# Ensure Conda is available    
if [ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
    echo "Conda not found at ${CONDA_BASE}. Please check Conda installation."
    exit 1
fi
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate PREDICT

# === Example for each module ===
# For the meaning of each argument, please see README.md "Arguments".

# --- FindKmers ---
# Basic usage (use default arguments)
python /path/to/PREDICT/PREDICT.py FindKmers

# Custom arguments
# python /path/to/PREDICT/PREDICT.py FindKmers --features gene --pThreshold 0.01 --alg RandomForest --up_stream 2000 --down_stream 800


# --- Kmer2Motif ---
# Basic usage (use default arguments)
# python /path/to/PREDICT/PREDICT.py Kmer2Motif

# Custom arguments
# python /path/to/PREDICT/PREDICT.py Kmer2Motif --TopKmers 0.05 --KeepTopMotifs 0.2 --ScoreCutOff 0.85


# --- RanKmers ---
# Basic usage (use default arguments)
# python /path/to/PREDICT/PREDICT.py RanKmers

# Custom arguments
# python /path/to/PREDICT/PREDICT.py RanKmers --features gene --alg RandomForest --up_stream 300 --down_stream 700


# --- TeamKmers ---
# Basic usage (use default arguments)
# python /path/to/PREDICT/PREDICT.py TeamKmers

# Custom arguments
# python /path/to/PREDICT/PREDICT.py TeamKmers --Motif False --TPCSThreshold 0.75 --TPTNCSThreshold 0.3


# --- ViewKmers --- NOT recommended to run this module on HPC, better to run it locally in terminal.
# Basic usage (use default arguments)
# streamlit run /path/to/PREDICT/PREDICT.py -- ViewKmers

# Custom arguments
# streamlit run /path/to/PREDICT/PREDICT.py -- ViewKmers --features gene --up_stream 1000 --down_stream 500



