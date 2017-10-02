# to run this script: 

system("make");

# C17
system("./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C17_mapped.blif ./iscas/C17out200.pl ./results/C17ED.txt > ./results/C17_450.txt");
print "C17 finished!\n";

# C432
system("./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C432_mapped.blif ./iscas/C432out1200w.pl ./results/C432ED.txt > ./results/C432_450.txt");
print "C432 finished!\n";

# C499
system("./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C499_mapped.blif ./iscas/C499out2500w.pl ./results/C499ED.txt > ./results/C499_450.txt");
print "C499 finished!\n";

# C880
system("./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C880_mapped.blif ./iscas/C880out1700w.pl ./results/C880ED.txt > ./results/C880_450.txt");
print "C880 finished!\n";

# C1355
system("./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C1355_mapped.blif ./iscas/C1355out2500w.pl ./results/C1355ED.txt > ./results/C1355_450.txt");
print "C1355 finished!\n";

# C1908
system("./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C1908_mapped.blif ./iscas/C1908out2100w.pl ./results/C1908ED.txt > ./results/C1908_450.txt");
print "C1908 finished!\n";

# C2670
system("./asloptimizer ./iscas/10nm_ASL.genlib ./iscas/C2670_mapped.blif ./iscas/C2670out2500.pl ./results/C2670ED.txt > ./results/C2670_450.txt");
print "C2670 finished!\n";

# C3540
system("./asloptimizer ../iscas/10nm_ASL.genlib ../iscas/C3540_mapped.blif ../iscas/C3540out3100.pl ./results/C3540ED.txt > ./results/C3540_450.txt");
print "C3540 finished!\n";

# C5315
system("./asloptimizer ../iscas/10nm_ASL.genlib ../iscas/C5315_mapped.blif ../iscas/C5315out3600.pl ./results/C5315ED.txt > ./results/C5315_450.txt");
print "C5315 finished!\n";

# C6288
system("./asloptimizer ../iscas/10nm_ASL.genlib ../iscas/C6288_mapped.blif ../iscas/C6288out5100.pl ./results/C6288ED.txt > ./results/C6288_450.txt");
print "C6288 finished!\n";

# C7552
system("./asloptimizer ../iscas/10nm_ASL.genlib ../iscas/C7552_mapped.blif ../iscas/C7552out4000.pl ./results/C7552ED.txt > ./results/C7552_450.txt");
print "C7552 finished!\n";








