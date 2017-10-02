#!/bin/bash
# A bash script to run the optimizer on all the ISCAS-85 benchmarks in ./iscas folder

echo "Executing script by Z. L. on $(hostname) at $(date)"

echo "Begin to compile asloptimizer: "
make
echo "Compilation success"

mkdir -p ./results

echo "Optimizing 17"
./asloptimizer ../iscas/10nm_ASL.genlib ../iscas/C17_mapped.blif ../iscas/C17out200.pl ./results/C17ED.txt > ./results/C17_450.txt
echo "C17 optimization finished!"


echo "Removing execution file and .o files: "

make clean

echo "Finish optimization on ISCAS benchmarks at $(date)" 
