# include "gate.h"

void GateDelayCalc(vector<gate*> &gate_vec, map<string, gate*> &gate_map)
{
	for (int i=0; i<gate_vec.size(); i++)
	{
		IntDelayCalc(gate_vec[i]); // calculate the intrinsic delay of each gate
		
		gate_vec[i]->ext_ipt_delay.resize(gate_vec[i]->ipt_ptrs.size()); // calculate the external delay from each input of this gate gate_vec[i]
		for (int j=0; j<gate_vec[i]->ipt_ptrs.size(); j++)
		{
			InterconnectDelayCalc(gate_vec[i]->ipt_ptrs[j], gate_vec[i], j); // calculate the delay from input gate to output gate
		}

		if (gate_vec[i]->gatetype == "PRI_INPUT")
		{
			gate_vec[i]->RTready = 1;
			// ReachTime of Primary input is initialized as 0
			gate_vec[i]->ReachTime = 0;
		}
		else
		{
			gate_vec[i]->RTready = 0;
			// ReachTime of each gate is initialized its intrinsic delay
			gate_vec[i]->ReachTime = gate_vec[i]->int_ipt_delay[0];
		}

	}
}

/* this version should only be used after initialization (means after the version GateDelayCalc without gate *sized_gate_ptr (sgate_ptr) is executed) */
void GateDelayCalc(vector<gate*> &gate_vec, map<string, gate*> &gate_map, gate *sgate_ptr)
{
//	cout<<"Internal delay for "<<sgate_ptr->name<<" is recalculated!"<<endl;
	IntDelayCalc(sgate_ptr); // only calculate the intrinsic delay of sized gate

	// gate_vec[i]->ext_ipt_delay.resize(gate_vec[i]->ipt_ptrs.size()); // no need to resize the delay storage vector // calculate the external delay from each input of this gate gate_vec[i]
	// only recalcualte the external input delay of that sized gate
	//	cout<<"The external input interconnect delay of "<<sgate_ptr->name<<" is recalculated!"<<endl;
	for (int i=0; i<sgate_ptr->ipt_ptrs.size(); i++)
	{
	//	cout<<"The external input interconnect delay from "<<gate_vec[i]->ipt_ptrs[j]->name<<" to "<<sgate_ptr->name<<" is recalculated!"<<endl;
		InterconnectDelayCalc(sgate_ptr->ipt_ptrs[i], sgate_ptr, i); // calculate the delay from input gate to output gate
	}

	gate *opt_sgate_ptr;
	for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
	{
		opt_sgate_ptr = sgate_ptr->opt_ptrs[i];
		int position = std::find(opt_sgate_ptr->ipt_ptrs.begin(), opt_sgate_ptr->ipt_ptrs.end(), sgate_ptr) - opt_sgate_ptr->ipt_ptrs.begin();
		InterconnectDelayCalc(sgate_ptr,opt_sgate_ptr, position);
	}

	for (int i=0; i<gate_vec.size(); i++)
	{
		// this part of code is commented if we need to perform incremental timing analysis, which needs the unimpacted gates still keep the RTready == 1
		if (gate_vec[i]->gatetype == "PRI_INPUT")
		{
			gate_vec[i]->RTready = 1;
			// ReachTime of Primary input is initialized as 0
			gate_vec[i]->ReachTime = 0;
		}
		else
		{
			gate_vec[i]->RTready = 0;
			// ReachTime of each gate is initialized its intrinsic delay
			gate_vec[i]->ReachTime = gate_vec[i]->int_ipt_delay[0];
		}
	}

}

/* this version further reduced the number of input parameters and should be only used in versions with incremental timing analysis */
void GateDelayCalc(gate *sgate_ptr)
{
//	cout<<"Internal delay for "<<sgate_ptr->name<<" is recalculated!"<<endl;
	IntDelayCalc(sgate_ptr); // only calculate the intrinsic delay of sized gate

	// gate_vec[i]->ext_ipt_delay.resize(gate_vec[i]->ipt_ptrs.size()); // no need to resize the delay storage vector // calculate the external delay from each input of this gate gate_vec[i]
	// only recalcualte the external input delay of that sized gate
	//	cout<<"The external input interconnect delay of "<<sgate_ptr->name<<" is recalculated!"<<endl;
	for (int i=0; i<sgate_ptr->ipt_ptrs.size(); i++)
	{
	//	cout<<"The external input interconnect delay from "<<gate_vec[i]->ipt_ptrs[j]->name<<" to "<<sgate_ptr->name<<" is recalculated!"<<endl;
		InterconnectDelayCalc(sgate_ptr->ipt_ptrs[i], sgate_ptr, i); // calculate the delay from input gate to output gate
	}

	gate *opt_sgate_ptr;
	for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
	{
		opt_sgate_ptr = sgate_ptr->opt_ptrs[i];
		int position = std::find(opt_sgate_ptr->ipt_ptrs.begin(), opt_sgate_ptr->ipt_ptrs.end(), sgate_ptr) - opt_sgate_ptr->ipt_ptrs.begin();
		InterconnectDelayCalc(sgate_ptr, opt_sgate_ptr, position);
	}

}

/* this version further uses a precalculation (lookup) table to find the Internal delay */
void GateDelayCalc(gate *sgate_ptr, map<pair<int, int>, double> &Delay3to1Table, map<pair<int, int>, double> &Delay5to1Table)
{
//	cout<<"Internal delay for "<<sgate_ptr->name<<" is recalculated!"<<endl;
//	IntDelayCalc(sgate_ptr); // only calculate the intrinsic delay of sized gate
	IntDelayCalcMatch(sgate_ptr, Delay3to1Table, Delay5to1Table);

	// gate_vec[i]->ext_ipt_delay.resize(gate_vec[i]->ipt_ptrs.size()); // no need to resize the delay storage vector // calculate the external delay from each input of this gate gate_vec[i]
	// only recalcualte the external input delay of that sized gate
	//	cout<<"The external input interconnect delay of "<<sgate_ptr->name<<" is recalculated!"<<endl;
	for (int i=0; i<sgate_ptr->ipt_ptrs.size(); i++)
	{
	//	cout<<"The external input interconnect delay from "<<gate_vec[i]->ipt_ptrs[j]->name<<" to "<<sgate_ptr->name<<" is recalculated!"<<endl;
		InterconnectDelayCalc(sgate_ptr->ipt_ptrs[i], sgate_ptr, i); // calculate the delay from input gate to output gate
	}

	gate *opt_sgate_ptr;
	for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
	{
		opt_sgate_ptr = sgate_ptr->opt_ptrs[i];
		int position = std::find(opt_sgate_ptr->ipt_ptrs.begin(), opt_sgate_ptr->ipt_ptrs.end(), sgate_ptr) - opt_sgate_ptr->ipt_ptrs.begin();
		InterconnectDelayCalc(sgate_ptr, opt_sgate_ptr, position);
	}

}

/* this version adds parameters to be extracted and used in recovery */
void GateDelayExtractNCalc(gate *sgate_ptr, double &TempInternalDelay, vector<double> &TempInputDelay, vector<double> &TempOutputDelay)
{
	if (sgate_ptr->int_ipt_delay.size() > 0)
	{
		TempInternalDelay = sgate_ptr->int_ipt_delay[0];
		IntDelayCalc(sgate_ptr); // only calculate the intrinsic delay of sized gate
	}

	// gate_vec[i]->ext_ipt_delay.resize(gate_vec[i]->ipt_ptrs.size()); // no need to resize the delay storage vector // calculate the external delay from each input of this gate gate_vec[i]
	// only recalcualte the external input delay of that sized gate
	//	cout<<"The external input interconnect delay of "<<sgate_ptr->name<<" is recalculated!"<<endl;
	for (int i=0; i<sgate_ptr->ipt_ptrs.size(); i++)
	{
	//	cout<<"The external input interconnect delay from "<<gate_vec[i]->ipt_ptrs[j]->name<<" to "<<sgate_ptr->name<<" is recalculated!"<<endl;
		TempInputDelay.push_back(sgate_ptr->ext_ipt_delay[i]);
		InterconnectDelayCalc(sgate_ptr->ipt_ptrs[i], sgate_ptr, i); // calculate the delay from input gate to output gate
	}

	gate *opt_sgate_ptr;
	for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
	{
		opt_sgate_ptr = sgate_ptr->opt_ptrs[i];
		int position = std::find(opt_sgate_ptr->ipt_ptrs.begin(), opt_sgate_ptr->ipt_ptrs.end(), sgate_ptr) - opt_sgate_ptr->ipt_ptrs.begin();
		TempOutputDelay.push_back(opt_sgate_ptr->ext_ipt_delay[position]);
		InterconnectDelayCalc(sgate_ptr, opt_sgate_ptr, position);
	}

}

/* this version adds parameters to be extracted and used in recovery, and use lookup table to find internal delay */
void GateDelayExtractNCalc(gate *sgate_ptr, double &TempInternalDelay, vector<double> &TempInputDelay, vector<double> &TempOutputDelay, map<pair<int, int>, double> &Delay3to1Table, map<pair<int, int>, double> &Delay5to1Table)
{
	if (sgate_ptr->int_ipt_delay.size() > 0)
	{
		TempInternalDelay = sgate_ptr->int_ipt_delay[0];
//		IntDelayCalc(sgate_ptr); // only calculate the intrinsic delay of sized gate
		IntDelayCalcMatch(sgate_ptr, Delay3to1Table, Delay5to1Table);
	//	if (sgate_ptr->name == "n9") cout<<sgate_ptr->int_ipt_delay[0]<<endl;
	}

	// gate_vec[i]->ext_ipt_delay.resize(gate_vec[i]->ipt_ptrs.size()); // no need to resize the delay storage vector // calculate the external delay from each input of this gate gate_vec[i]
	// only recalcualte the external input delay of that sized gate
	//	cout<<"The external input interconnect delay of "<<sgate_ptr->name<<" is recalculated!"<<endl;
	for (int i=0; i<sgate_ptr->ipt_ptrs.size(); i++)
	{
	//	cout<<"The external input interconnect delay from "<<gate_vec[i]->ipt_ptrs[j]->name<<" to "<<sgate_ptr->name<<" is recalculated!"<<endl;
		TempInputDelay.push_back(sgate_ptr->ext_ipt_delay[i]);
		InterconnectDelayCalc(sgate_ptr->ipt_ptrs[i], sgate_ptr, i); // calculate the delay from input gate to output gate
	}

	gate *opt_sgate_ptr;
	for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
	{
		opt_sgate_ptr = sgate_ptr->opt_ptrs[i];
		int position = std::find(opt_sgate_ptr->ipt_ptrs.begin(), opt_sgate_ptr->ipt_ptrs.end(), sgate_ptr) - opt_sgate_ptr->ipt_ptrs.begin();
		TempOutputDelay.push_back(opt_sgate_ptr->ext_ipt_delay[position]);
		InterconnectDelayCalc(sgate_ptr, opt_sgate_ptr, position);
	}

}

/* this version further select the fan-in or fan-out delay to be recalculated based on the sized magnet is input or output */
void GateDelayExtractNCalc(gate *sgate_ptr, string ioMagnetSelect, double &TempInternalDelay, vector<double> &TempInputDelay, vector<double> &TempOutputDelay, map<pair<int, int>, double> &Delay3to1Table, map<pair<int, int>, double> &Delay5to1Table)
{
	if (sgate_ptr->int_ipt_delay.size() > 0)
	{
		TempInternalDelay = sgate_ptr->int_ipt_delay[0];
//		IntDelayCalc(sgate_ptr); // only calculate the intrinsic delay of sized gate
		IntDelayCalcMatch(sgate_ptr, Delay3to1Table, Delay5to1Table);
	//	if (sgate_ptr->name == "n9") cout<<sgate_ptr->int_ipt_delay[0]<<endl;
	}

	// gate_vec[i]->ext_ipt_delay.resize(gate_vec[i]->ipt_ptrs.size()); // no need to resize the delay storage vector // calculate the external delay from each input of this gate gate_vec[i]
	// only recalcualte the external input delay of that sized gate
	//	cout<<"The external input interconnect delay of "<<sgate_ptr->name<<" is recalculated!"<<endl;

	if ( (ioMagnetSelect != "input") && (ioMagnetSelect != "output")) {	cout<<"Error: Undefined ioMagnetSelect in fucntion GateDelayExtractNCalc!!!"<<endl; assert(false);	}

	for (int i=0; i<sgate_ptr->ipt_ptrs.size(); i++)
	{
	//	cout<<"The external input interconnect delay from "<<gate_vec[i]->ipt_ptrs[j]->name<<" to "<<sgate_ptr->name<<" is recalculated!"<<endl;
		TempInputDelay.push_back(sgate_ptr->ext_ipt_delay[i]);
		if ((ioMagnetSelect == "input") || (std::find(gatetype1to1list.begin(), gatetype1to1list.end(), sgate_ptr->gatetype) != gatetype1to1list.end()))
		{
			InterconnectDelayCalc(sgate_ptr->ipt_ptrs[i], sgate_ptr, i); // calculate the delay from input gate to output gate
		}		
	}

	gate *opt_sgate_ptr;
	for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
	{
		opt_sgate_ptr = sgate_ptr->opt_ptrs[i];
		int position = std::find(opt_sgate_ptr->ipt_ptrs.begin(), opt_sgate_ptr->ipt_ptrs.end(), sgate_ptr) - opt_sgate_ptr->ipt_ptrs.begin();
		TempOutputDelay.push_back(opt_sgate_ptr->ext_ipt_delay[position]);
		if (ioMagnetSelect == "output")
		{
			InterconnectDelayCalc(sgate_ptr, opt_sgate_ptr, position);
		}
	}
}

/* this version adds parameters to be extracted and used in recovery, and extract the related sized delays as well */
void GateDelayExtractNCalc(gate *sgate_ptr, double &TempInternalDelay, vector<double> &TempInputDelay, vector<double> &TempOutputDelay, double &SizedInternalDelay, vector<double> &SizedInputDelay, vector<double> &SizedOutputDelay)
{
	if (sgate_ptr->int_ipt_delay.size() > 0)
	{
		TempInternalDelay = sgate_ptr->int_ipt_delay[0];
		IntDelayCalc(sgate_ptr); // only calculate the intrinsic delay of sized gate
		SizedInternalDelay = sgate_ptr->int_ipt_delay[0];
	}

	// gate_vec[i]->ext_ipt_delay.resize(gate_vec[i]->ipt_ptrs.size()); // no need to resize the delay storage vector // calculate the external delay from each input of this gate gate_vec[i]
	// only recalcualte the external input delay of that sized gate
	//	cout<<"The external input interconnect delay of "<<sgate_ptr->name<<" is recalculated!"<<endl;
	for (int i=0; i<sgate_ptr->ipt_ptrs.size(); i++)
	{
	//	cout<<"The external input interconnect delay from "<<gate_vec[i]->ipt_ptrs[j]->name<<" to "<<sgate_ptr->name<<" is recalculated!"<<endl;
		TempInputDelay.push_back(sgate_ptr->ext_ipt_delay[i]);
		InterconnectDelayCalc(sgate_ptr->ipt_ptrs[i], sgate_ptr, i); // calculate the delay from input gate to output gate
		SizedInputDelay.push_back(sgate_ptr->ext_ipt_delay[i]);
	}

	gate *opt_sgate_ptr;
	for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
	{
		opt_sgate_ptr = sgate_ptr->opt_ptrs[i];
		int position = std::find(opt_sgate_ptr->ipt_ptrs.begin(), opt_sgate_ptr->ipt_ptrs.end(), sgate_ptr) - opt_sgate_ptr->ipt_ptrs.begin();
		TempOutputDelay.push_back(opt_sgate_ptr->ext_ipt_delay[position]);
		InterconnectDelayCalc(sgate_ptr, opt_sgate_ptr, position);
		SizedOutputDelay.push_back(opt_sgate_ptr->ext_ipt_delay[position]);
	}

}

/* this function works together with GateDelayExtractNCalc to recover the orginial delays related to the sized gate */
void GateDelayRecover(gate *sgate_ptr, double &InternalDelay, vector<double> &InputDelay, vector<double> &OutputDelay)
{
	if (sgate_ptr->int_ipt_delay.size() > 0)
	{
		fill(sgate_ptr->int_ipt_delay.begin(), sgate_ptr->int_ipt_delay.end(), InternalDelay);
	}

	for (int i=0; i<sgate_ptr->ipt_ptrs.size(); i++)
	{
		sgate_ptr->ext_ipt_delay[i] = InputDelay[i];
	}

	gate *opt_sgate_ptr;
	for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
	{
		opt_sgate_ptr = sgate_ptr->opt_ptrs[i];
		int position = std::find(opt_sgate_ptr->ipt_ptrs.begin(), opt_sgate_ptr->ipt_ptrs.end(), sgate_ptr) - opt_sgate_ptr->ipt_ptrs.begin();
		opt_sgate_ptr->ext_ipt_delay[position] = OutputDelay[i];
	}

}

/* Extract the delays related to the sized gate on the critical path */
double ExtractSizedGateRelatedDelays(gate *sgate_ptr, vector<gate*> &CriticalPathGates)
{
	double SizedGateRelatedDelays = 0;
	if (sgate_ptr->gatetype == "PRI_INPUT")
	{
		// for the delay change of the fan-out gates of this sized gate on the critical path (for primary input gate, this is the only delay it influences)
		for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
		{
			gate *ogate_ptr;		
			if ( std::find(CriticalPathGates.begin(), CriticalPathGates.end(), sgate_ptr->opt_ptrs[i]) != CriticalPathGates.end() )
			{
				ogate_ptr = sgate_ptr->opt_ptrs[i];
				int position = std::find(ogate_ptr->ipt_ptrs.begin(), ogate_ptr->ipt_ptrs.end(), sgate_ptr) - ogate_ptr->ipt_ptrs.begin();
				SizedGateRelatedDelays += ogate_ptr->ext_ipt_delay[position];
			}
		}
	}
	else if ( std::find(gatetype1to1list.begin(), gatetype1to1list.end(), sgate_ptr->gatetype) != gatetype1to1list.end() )
	{
		// for the delay change of the fan-in gates of this sized gate on the critical path
		SizedGateRelatedDelays += sgate_ptr->ext_ipt_delay[0];

		// for the delay change of the fan-out gates of this sized gate on the critical path
		gate *ogate_ptr;
		for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
		{
			if ( std::find(CriticalPathGates.begin(), CriticalPathGates.end(), sgate_ptr->opt_ptrs[i]) != CriticalPathGates.end() )
			{
				ogate_ptr = sgate_ptr->opt_ptrs[i];
				int position = std::find(ogate_ptr->ipt_ptrs.begin(), ogate_ptr->ipt_ptrs.end(), sgate_ptr) - ogate_ptr->ipt_ptrs.begin();
				SizedGateRelatedDelays += ogate_ptr->ext_ipt_delay[position];
			}
		}		
	}
	else if ( ( std::find(gatetype3to1list.begin(), gatetype3to1list.end(), sgate_ptr->gatetype) != gatetype3to1list.end() ) || ( std::find(gatetype5to1list.begin(), gatetype5to1list.end(), sgate_ptr->gatetype) != gatetype5to1list.end() ) )
	{
		// for the delay change of the fan-in gates and the internal delay of this sized gate on the critical path
		gate *igate_ptr;
		for (int i=0; i<sgate_ptr->ipt_ptrs.size(); i++)
		{
			if ( std::find(CriticalPathGates.begin(), CriticalPathGates.end(), sgate_ptr->ipt_ptrs[i]) != CriticalPathGates.end() )
			{	
				igate_ptr = sgate_ptr->ipt_ptrs[i];
				int position = std::find(sgate_ptr->ipt_ptrs.begin(), sgate_ptr->ipt_ptrs.end(), igate_ptr) - sgate_ptr->ipt_ptrs.begin();
				SizedGateRelatedDelays += sgate_ptr->int_ipt_delay[position];
				SizedGateRelatedDelays += sgate_ptr->ext_ipt_delay[position];
			}
		}

		// for the delay change of the fan-out gates of this sized gate on the critical path
		gate *ogate_ptr;
		for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
		{
			if ( std::find(CriticalPathGates.begin(), CriticalPathGates.end(), sgate_ptr->opt_ptrs[i]) != CriticalPathGates.end() )
			{
				ogate_ptr = sgate_ptr->opt_ptrs[i];
				int position = std::find(ogate_ptr->ipt_ptrs.begin(), ogate_ptr->ipt_ptrs.end(), sgate_ptr) - ogate_ptr->ipt_ptrs.begin();
				SizedGateRelatedDelays += ogate_ptr->ext_ipt_delay[position];
			}
		}


	}
	else
	{	
		cout<<"Error: Undefined gatetype in function ExtractSizedGateRelatedDelays!!!"<<endl;
		assert(false);
	}

	return SizedGateRelatedDelays;
}


