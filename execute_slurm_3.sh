#!/bin/bash
#SBATCH --job-name=n_bod
#SBATCH --partition=Centaurus
#SBATCH --time=03:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256G
./run_3.sh
