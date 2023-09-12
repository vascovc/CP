clc
clearvras
close all

% Este script MATLAB permite importar o 
% ficheiro bina'rio produzido pelos programas em C
% e produzir um plot. Pode ser muito util para
% ver se as condicoes fronteira estao corretas,
% para ve se nao ha' erros nas fronteiras dos
% subdominios, etc.

M = 10;

% fileID = fopen('output_v1.bin');
fileID = fopen('output_v2.bin');

a = fread(fileID, M, 'double');
b = fread(fileID, M, 'int');
fclose(fileID);

