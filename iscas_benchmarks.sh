#!/bin/bash
# A bash script to run the optimizer on all the ISCAS-85 benchmarks in ./iscas folder

echo "Executing script by Z. L. on $(hostname) at $(date)"
echo
echo "Begin to compile asloptimizer:"
make
echo "Compilation success!"
echo
mkdir -p ./results

echo "Optimizing C17"
./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C17_mapped.blif ./iscas/C17out200.pl ./results/C17ED.txt > ./results/C17_450.txt
echo "C17 optimization finished!"
echo

echo "Optimizing C432"
./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C432_mapped.blif ./iscas/C432out1200w.pl ./results/C432ED.txt > ./results/C432_450.txt
echo "C432 optimization finished!"
echo

echo "Optimizing C499"
./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C499_mapped.blif ./iscas/C499out2500w.pl ./results/C499ED.txt > ./results/C499_450.txt
echo "C499 optimization finished!"
echo

echo "Optimizing C880"
./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C880_mapped.blif ./iscas/C880out1700w.pl ./results/C880ED.txt > ./results/C880_450.txt
echo "C880 optimization finished!"
echo

echo "Optimizing C1355"
./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C1355_mapped.blif ./iscas/C1355out2500w.pl ./results/C1355ED.txt > ./results/C1355_450.txt
echo "C1355 optimization finished!"
echo

echo "Optimizing C1908"
./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C1908_mapped.blif ./iscas/C1908out2100w.pl ./results/C1908ED.txt > ./results/C1908_450.txt
echo "C1908 optimization finished!"
echo

echo "Optimizing C2670"
./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C2670_mapped.blif ./iscas/C2670out2500.pl ./results/C2670ED.txt > ./results/C2670_450.txt
echo "C2670 optimization finished!"
echo

echo "Optimizing C3540"
./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C3540_mapped.blif ./iscas/C3540out3100.pl ./results/C3540ED.txt > ./results/C3540_450.txt
echo "C3540 optimization finished!"
echo

echo "Optimizing C5315"
./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C5315_mapped.blif ./iscas/C5315out3600.pl ./results/C5315ED.txt > ./results/C5315_450.txt
echo "C5315 optimization finished!"
echo

echo "Optimizing C6288"
./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C6288_mapped.blif ./iscas/C6288out5100.pl ./results/C6288ED.txt > ./results/C6288_450.txt
echo "C6288 optimization finished!"
echo

echo "Optimizing C7552"
./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C7552_mapped.blif ./iscas/C7552out4000.pl ./results/C7552ED.txt > ./results/C7552_450.txt
echo "C7552 optimization finished!"
echo


echo "Removing execution file and .o files:"
make clean
echo
echo "Finish optimization on ISCAS benchmarks at $(date)" 
