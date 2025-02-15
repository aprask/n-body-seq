#!/bin/bash
make clean
make
if [ $? -ne 0 ]; then
    echo "Build failed. Exiting."
    exit 1
fi 
./nbod 100 1 10000 500
