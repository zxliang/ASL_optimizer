# ASL_optimizer
An optimizer for perform TILOS-like algorithm for energy/delay trade-offs in ASL circuits
==========================

Author: Zhaoxin Liang @ University of Minnesota - Twin Cities, VLSI - EDA group
Description: This repository contains C++ souce codes for the ASL geometric optimizer
Version: 1.0
License: by Zhaoxin Liang

Note 1: Created for ICCAD 2015 submission
Note 2: Version 1.0 uses Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page for matrix computation in MNA.cpp.
Note 3: test sync.


Note 4: v6 does not include incremental timing analysis, it is included in v5
Note 5: v7 removed the STA performance in every TILOS iterations.
Note 6: v8 added incremental timing analysis back
Note 7: The "same delay reduction flag" is removed in v9, yet the STA is recovered for the evaluation after each adaption.
Note 8: In version 10, a bug in incremental timing analysis by adding the "gate_ptr->opt_ptrs[i]->RTready != 0" in "ReachTimeReadyReset" function.
Note 9: In version 11, removed the recalculation and use recovery from stored value to save runtime in each iteration.
Note 10: In version 12, the adaption should also reuse previous calculation results.
	(it seems that the double precision matching issue with map<pair<double, double>, double> Delay3to1Table, meaning the value you calculated for double in precalculation and later matching function cannot match well. trying to use pair<int, int> instead)
Note 11: In version 13, recalculation of delays related to the sized gates is adjusted, so some unnecessary recalculations are removed, for example, if the output magnet is sized, only the internal delay and its fan-out interconnect delay should be recalculated, no need to recaclulate the input interconnect delays of this gate.
	(Finished. the most time-consuming part is the calculation of the gate/interconnect delays, so precalculation/reuse/avoiding unnecessary calculation are good ways to reduce runtime.)

Note 12 (Aug. 11, 2015): This program is used to generate data in JxCDC manuscript in Aug. 8, 2015. Three sets of technology parameters were used with the [double SpacingDis (from mylib.cpp)] parameter set to the same value as the corresponding spin diffusion lengths in the gate.h. For fast simulation speed, calculation and output of Power and Area is disabled by commenting function [vector<double> PowerArea = PowerAreaCalc(gate_vec0);] in row 653 of mylib.cpp and turn off the output line a few rows later [outfile<<CURRENT_ITERATION+1<<" "<<MaxDelay0<<" "<<PowerArea[0]<<" "<<PowerArea[1]<<endl;]. 


Note 13 (Oct. 16, 2015): New models considering the contacts and grounds were added to the delay calculation part.
