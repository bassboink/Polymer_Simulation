#!/bin/bash

cd ~/Desktop/TEK8/homo

mpiexec -np 4 lammps-daily < imdinputH.txt

