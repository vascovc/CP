#!/bin/bash

echo -e "\nbuilding.pgm"
for i in {1..10}
do
./fastDetectorCuda_better -i building.pgm
./testDiffs referenceCuda.pgm resultCuda.pgm
done

echo -e "\nchessBig.pgm"
for i in {1..10}
do
./fastDetectorCuda_better -i chessBig.pgm
./testDiffs referenceCuda.pgm resultCuda.pgm
done

echo -e "\nchessRotate1.pgm"
for i in {1..10}
do
./fastDetectorCuda_better -i chessRotate1.pgm
./testDiffs referenceCuda.pgm resultCuda.pgm
done

echo -e "\nhouse.pgm"
for i in {1..10}
do
./fastDetectorCuda_better -i house.pgm
./testDiffs referenceCuda.pgm resultCuda.pgm
done

echo -e "\nsquares.pgm"
for i in {1..10}
do
./fastDetectorCuda_better -i squares.pgm
./testDiffs referenceCuda.pgm resultCuda.pgm
done
