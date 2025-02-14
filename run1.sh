#!/bin/bash
make clean
make
if [ $? -ne 0 ]; then
    echo "Build failed. Exiting."
    exit 1
fi 
./nbod solar 200 5000000 10000
