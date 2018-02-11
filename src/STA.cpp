#include "gate.h"

/* this .cpp inlcudes functions to perform STA and energy calculation */

/* Incremental Timing Analysis is used in optimization to reduce run time */
double IncrementalTimingAnalysis(vector<gate*> &gate_vec, map<string, gate*> &gate_map, vector<gate*> &CriticalPathGates, gate *sgate_ptr)
{
	ReachTimeReadyReset(sgate_ptr);
//	PrintGateReachTime(gate_vec);
	
	double MaxDelay = 0;
	for (int k=0; k<gate_vec.size(); k++)
	{
		gate_vec[k]->ReachTime = FindReachTime(gate_vec[k]);
		if ( (gate_vec[k]->RTready == 1) && (MaxDelay<gate_vec[k]->ReachTime) )
		{
			MaxDelay = gate_vec[k]->ReachTime;
		}
	}

	for (int l=0; l<gate_vec.size(); l++)
	{
		if (gate_vec[l]->ReachTime == MaxDelay)
		{
			CriticalPathGates.push_back(gate_vec[l]);
			FindCriticalPathGates(gate_vec[l], CriticalPathGates);
		}
	}

//	PrintGateReachTime(gate_vec);
	return MaxDelay;
}

/* reset the RTready value starting from the changed gate, for incremental timing analysis use */
void ReachTimeReadyReset(gate *gate_ptr)
{
	// this function requires all the RTready flag in gates all be set to 1 already, which is enabled in GateDelayCalc 
	// function with sizing gate pointer input version
//	cout<<gate_ptr->name<<endl;
	// if this gate is not primary input, set RTready to 0
	if (gate_ptr->gatetype != "PRI_INPUT")
	{
		gate_ptr->RTready = 0;
	}
	
	// set RTready value of the fan-outs of this gate to 0 recursively
	for (int i=0; i<gate_ptr->opt_ptrs.size(); i++)
	{
		if (gate_ptr->opt_ptrs[i]->RTready != 0)
		{
			ReachTimeReadyReset(gate_ptr->opt_ptrs[i]);
		}		
	}
}

double StaticTimingAnalysis(vector<gate*> &gate_vec, map<string, gate*> &gate_map, vector<gate*> &CriticalPathGates)
{
	double MaxDelay = 0;
	for (int k=0; k<gate_vec.size(); k++)
	{
		gate_vec[k]->ReachTime = FindReachTime(gate_vec[k]);
		if ( (gate_vec[k]->RTready == 1) && (MaxDelay<gate_vec[k]->ReachTime) )
		{
			MaxDelay = gate_vec[k]->ReachTime;
		}
	}

	for (int l=0; l<gate_vec.size(); l++)
	{
		if (gate_vec[l]->ReachTime == MaxDelay)
		{
			CriticalPathGates.push_back(gate_vec[l]);
			FindCriticalPathGates(gate_vec[l], CriticalPathGates);
		}
	}

	return MaxDelay;
}

double FindReachTime(gate *gate_ptr)
{
	double MaxReachTime = 0;

	if (gate_ptr->RTready == 1)
	{
		return gate_ptr->ReachTime;
	}
	else if (gate_ptr->RTready == 0)
	{
		gate_ptr->InputReachTime.resize(gate_ptr->ipt_ptrs.size());

		for (int i=0; i<gate_ptr->InputReachTime.size(); i++)
		{
			gate_ptr->InputReachTime[i] = FindReachTime(gate_ptr->ipt_ptrs[i]) + gate_ptr->int_ipt_delay[i] + gate_ptr->ext_ipt_delay[i];	
		}

		MaxReachTime = *std::max_element(gate_ptr->InputReachTime.begin(), gate_ptr->InputReachTime.end());
		gate_ptr->ReachTime = MaxReachTime;
		gate_ptr->RTready = 1;

	}

	return MaxReachTime;
}

void FindCriticalPathGates(gate *gate_ptr, vector<gate*> &CriticalPathGates)
{
	double MaxInputDelay = *std::max_element(gate_ptr->InputReachTime.begin(), gate_ptr->InputReachTime.end());
	
	for (int m=0; m<gate_ptr->ipt_ptrs.size(); m++)
	{
		// if this input has the max InputReachTime AND this input gate is not a primary input (PRI_INPUT) since we do not size primary
		// input this time AND this input is not added to the CriticalPathGates yet, THEN we add it into the list */
		if ( (gate_ptr->InputReachTime[m] == MaxInputDelay) && (std::find(CriticalPathGates.begin(), CriticalPathGates.end(), gate_ptr->ipt_ptrs[m]) == CriticalPathGates.end() ) )
		{
			if (gate_ptr->ipt_ptrs[m]->gatetype != "PRI_INPUT")
			{
				CriticalPathGates.push_back(gate_ptr->ipt_ptrs[m]);
				FindCriticalPathGates(gate_ptr->ipt_ptrs[m], CriticalPathGates);
			}
			else if (gate_ptr->ipt_ptrs[m]->gatetype == "PRI_INPUT")
			{
				CriticalPathGates.push_back(gate_ptr->ipt_ptrs[m]);
			}
		}
	}
}


