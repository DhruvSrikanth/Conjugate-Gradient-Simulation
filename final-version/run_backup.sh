#!/bin/bash

# module load openmpi

mpic++ conjugate_gradient.cpp -o conjugate_gradient -O3 -ffast-math -march=native -mtune=native -lmpi

mpiexec --bind-to none --report-bindings conjugate_gradient 2560 parallel

rm conjugate_gradient

python3 visualize.py