#!/bin/bash

make

# correr individualmente cada um dos outros para o num de threads que se queiram com 
./run_cuda.sh > cuda.txt
./run_cuda_m.sh > cuda_m.txt
./run_cuda_better.sh > cuda_better.txt
./run_cuda_m_better.sh > cuda_m_better.txt

grep "Host processing time: " cuda*.txt > results_cuda_host.txt
grep "Device processing time: " cuda*.txt > results_cuda_device.txt