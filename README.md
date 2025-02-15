# Sequential N-Body Simulation in C++

## Introduction
The goal of this program is to measure the performance of a sequential N-Body simulation on Centaurus.

---

## Table of Contents
1. [Running](#running)
2. [Requirements](#requirements)
3. [Program](#program)
4. [Benchmark](#benchmark)
5. [References](#references)

---

## Running

- To enter project use the `cd <path>` command
- You can run the project by running `Make` followed by the name of the executable `./nbod <tsv file or N> <delta t> <time steps> <dump rate>`
- You can similarly run the program via the shell scripts. Example `./run.sh`. Please make sure to set permissions for exeuction via `chmod +x`
- In order to plot the file, you need to create a virtual environment, assuming Python 3 is installed. `python -m venv venv`
- Then install matplotlib via `pip install matplotlib` or `pip install -r requirements.txt`
- Then you can run the python script `python plot.py *.tsv data/*.pdf`
---

## Requirements
- g++ compiler
- Python 3
- matplotlib
- Centaurus

---

## Program
**Explanation:**

- Pairwise Calculation: I utilized a pair-wise calculation in the set of N bodies. This yields a time complexity of O(n^2) unlike Barnes-Hut which yields O(n log(n))
- Design: My original program utilized manual memory management w/ malloc and structs. I found this to become laborious. I went ahead and re-learned OOP in C++, as I don't use it very often. Vectors and objects made this program much more modular and easy to read. I tried to make this as memory efficient as possible, by passing by reference or pointer when possible.

---

## Benchmark
- 

---


## References
1. [Force Component](https://physics.stackexchange.com/questions/17285/split-gravitational-force-into-x-y-and-z-componenets)
2. [C++ Reference](https://cplusplus.com/reference/)

---

**Author:** Andrew T. Praskala 
**Date:** Februrary 15, 2025 
