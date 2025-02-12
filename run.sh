#!/bin/bash
make clean
make
if [ $? -ne 0 ]; then
    echo "Build failed. Exiting."
    exit 1
fi 
./nbod solar 200 5000000 5
python3 python3 plot.py results.tsv solar_output.pdf
./nbod 100 1 10000 5
python3 python3 plot.py results.tsv rand_output_1.pdf
./nbod 1000 1 10000 5
python3 python3 plot.py results.tsv rand_output_1.pdf