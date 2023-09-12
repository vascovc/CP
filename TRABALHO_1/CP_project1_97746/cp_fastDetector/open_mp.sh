#!/bin/bash

make -f Makefile_OpenMP

# correr individualmente cada um dos outros para o num de threads que se queiram com 
# ./run_open_mp.sh > open_mp_i_th.txt
# ./run_open_mp_m.sh > open_mp_m_i_th.txt

grep "Host processing time: " open_mp*.txt > results_open_mp_host.txt
grep "OpenMP processing time: " open_mp*.txt > results_open_mp_open.txt