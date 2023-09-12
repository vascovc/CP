close all
clear all
clc

nx = 100;
N = nx*nx;
ny = nx;

fileID = fopen('results_a.bin');

array_MPI = fread(fileID, [ny nx],'double');
fclose(fileID);

load("vnewMat.mat")

MSE = 1/N * sum((array_MPI-Vnew).^2,'all')