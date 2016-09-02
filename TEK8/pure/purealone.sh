#!/bin/bash

cd ~/Desktop/TEK8/pure

mpiexec -np 4 lammps-daily < imdinputpure.txt

