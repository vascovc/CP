#!/bin/bash

make
mpiexec -n 4 ./original
mpiexec -n 4 ./alinea_a
mpiexec -n 4 ./alinea_b
mpiexec -n 4 ./alinea_c
mpiexec -n 4 ./alinea_d