#!/bin/bash

cd ~/Desktop/TEK8/diblock

mpiexec -np 4 lammps-daily < imdinputdiblock.txt

