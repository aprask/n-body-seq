#!/bin/bash
#SBATCH --job-name=n_bod
#SBATCH --partition=Centaurus
#SBATCH --time=03:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128G
./run1.sh
./run2.sh
./run3.sh
