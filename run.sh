#!/bin/bash
make clean
make
if [ $? -ne 0 ]; then
    echo "Build failed. Exiting."
    exit 1
fi 
./nbod solar 200 5000000 5
./nbod 100 1 10000 5
./nbod 1000 1 10000 5