#include "gate.h"

double MaxDelay0;
double EnergyConsumption0;
double Power0;
double Area0;
vector<gate*> gate_vec0;
map<string, gate*> gate_map0;
vector<gate*> CriticalPathGates0;

map<pair<int, int>, double> Delay3to1TBL, Delay5to1TBL;

void Read_Blif_N_Pl(ifstream &bliff, ifstream &plf)
{

/* First section of this file read blif file and store the information to gate_vec0 and gate_map0 */	
	string line;
	int FLAG = 0;	// use this indicator to deal with multiple lines of input (later output maybe)
	while (getline(bliff,line)) 
	{
		if (line.length()<=1) continue;
		if (line[0]=='#') continue;
//		cout<<line<<endl;
		
		istringstream sword(line);
		string temp;
		while (sword>>temp) 
		{
			if (temp==".model") break;
			else if ( (temp==".inputs")||(FLAG==1) ) 
			{	
			//	cout<<WordCount(line)<<endl;
				if (FLAG != 1) 
				{	for (int i=0; i<WordCount(line)-1; i++) 
					{	sword>>temp;
						if (temp != "\\") 
						{
						//	cout<<temp<<endl;		
							gate *temp_gate_ptr = new gate;
							temp_gate_ptr->setGateInfo(temp,"PRI_INPUT",60,100);
							gate_vec0.push_back(temp_gate_ptr);
							gate_map0[temp] = temp_gate_ptr;
						
							FLAG = 0;
						}
						else if (temp == "\\") 
						{	
							FLAG = 1;
						}
					}
				}
				else if (FLAG == 1) 
				{	
					do{
						if (temp !="\\") 
						{
							gate *temp_gate_ptr = new gate;
							temp_gate_ptr->setGateInfo(temp,"PRI_INPUT",60,100);
							gate_vec0.push_back(temp_gate_ptr);
							gate_map0[temp] = temp_gate_ptr;
						
							FLAG = 0;
						}
						else if (temp == "\\") {
							FLAG = 1;
						}
					} while (sword>>temp);	

				}
			}
			else if (temp==".outputs") break;
			else if (temp==".end") break;
			else if (temp==".gate") 
			{ 
				gate *temp_gate_ptr = new gate;
				sword>>temp;
				string temptype=temp;
//				cout<<temptype<<" "<<WordCount(line)<<endl;
				for(int i=0; i<WordCount(line)-3; i++) 
				{
					sword>>temp;
					temp_gate_ptr->ipt_nodes.push_back(temp.substr(temp.find_first_of("=")+1));
//					cout<<temp<<" ";
				}
				sword>>temp;
//				cout<<temp<<endl;
				temp_gate_ptr->setGateInfo(temp.substr(temp.find_first_of("=")+1), temptype, 60, 100);

				gate_vec0.push_back(temp_gate_ptr);
				gate_map0[temp.substr(temp.find_first_of("=")+1)] = temp_gate_ptr;
			}
		}
	}

/* this part go through the whole pl find to find placement location information for each gate; 
at the same time, for each gate, the map between gate name and gate address is used to attach the input
gate address to the gate input list */

	FindLocation(plf, gate_map0);

/* this part for each gate, calculate the distances from one gate's input to this gate */

	DistanceCalc(gate_vec0, gate_map0);	

}


// attach the location information to each gate and find the gate address for each input node of that gate
void FindLocation(ifstream &plf, map<string, gate*> &gate_map)	
{
//	cout<<"finding nodes coordinates!"<<endl;
	string plline;
	while (getline(plf,plline)) 
	{
	//	cout<<plline<<endl;
		if (plline.length()<=1) continue;
		if (plline[0]=='#') continue;
		if (plline.substr(0,4)=="UCLA") continue;
	//	icout<<plline<<endl;
					
		istringstream snode(plline);
		string temppl;
		while (snode>>temppl) 
		{
			if ( gate_map.find(temppl) != gate_map.end() )
			{	
			//	cout<<temppl<<" gate found! with location x: ";
				gate *gate_ptr = gate_map.find(temppl)->second;
								
				snode>>temppl;
				stringstream(temppl) >> gate_ptr->optn_x;
				snode>>temppl;
				stringstream(temppl) >> gate_ptr->optn_y;

				gate_ptr->ipt_ptrs.resize(gate_ptr->ipt_nodes.size());
				for (int i=0; i<gate_ptr->ipt_nodes.size(); i++)
				{	
					gate *ipt_gate_ptr = new gate;
					ipt_gate_ptr = gate_map.find(gate_ptr->ipt_nodes[i])->second;
					gate_ptr->ipt_ptrs[i] = ipt_gate_ptr;
				}
		
	
				break;
			}
			else if ( gate_map.find(temppl) == gate_map.end() )
			{
			//	cout<<temppl<<" gate not found! Could be primary output gate!"<<endl;
				break;
			}
		
		}
	}
//	plf.clear();	// reading to the end of the file, eof flag needs to be cleared
//	plf.seekg(0, ios::beg);
}

// calculate the distance from each input as well as gate dimension initialization 
void DistanceCalc(vector<gate*> &gate_vec, map<string, gate*> &gate_map)
{
	for (int i=0; i<gate_vec.size(); i++)
	{	
		gate_vec[i]->ipt_dis.resize(gate_vec[i]->ipt_nodes.size());
		for (int j=0; j<gate_vec[i]->ipt_nodes.size(); j++)
		{
			gate_vec[i]->ipt_dis[j] = abs(gate_vec[i]->optn_x - (gate_vec[i]->ipt_ptrs[j])->optn_x) + abs(gate_vec[i]->optn_y - (gate_vec[i]->ipt_ptrs[j])->optn_y);
		}
		
		// gate dimension initialization
		gate_vec[i]->DimInit();
	}	
}


void Initialization(ofstream &outfile) 
{
	double SpacingDis = lambda_N*1e9; // 450/1000/2000 here use [nm] 
	BufferInsertion(gate_vec0, gate_map0, SpacingDis);
	AddOptPtrs(gate_vec0, gate_map0);

	GateDelayCalc(gate_vec0, gate_map0);
	MaxDelay0 = StaticTimingAnalysis(gate_vec0, gate_map0, CriticalPathGates0);
	EnergyConsumption0 = EnergyCalc(MaxDelay0, gate_vec0);
	Power0 = PowerCalc(gate_vec0);
	Area0 = AreaCalc(gate_vec0);

// this part do test print out the vecotor containing gate pointers 

	PrintGatePtrVec(gate_vec0);

	cout<<"MaxDelay: "<<MaxDelay0<<"ns"<<endl;
	cout<<"EnergyConsumption: "<<EnergyConsumption0<<"J"<<endl;
	cout<<"Power: "<<Power0<<"W"<<endl;
	cout<<"Area: "<<Area0<<"nm^2"<<endl;
	cout<<"Gates on Critical Path: ";
	
	for (int i=0; i<CriticalPathGates0.size(); i++)
	{
		cout<<CriticalPathGates0[i]->name<<" ";
	}
	cout<<endl;

	if (outfile.is_open() == 1)
	{
		outfile<<"0 "<<MaxDelay0<<" "<<Power0<<" "<<Area0<<endl;
	}

	InternalDelayPreCalc(geofctr, Delay3to1TBL, Delay5to1TBL, MagnetLengthUB);

}


/* This function inserted buffer (magnet) every SpacingDis (or closet to this distance) */
void BufferInsertion(vector<gate*> &gate_vec, map<string,gate*> &gate_map, double SpacingDis)
{
	for (int i=0; i<gate_vec.size(); i++)
	{
		for (int j=0; j<gate_vec[i]->ipt_nodes.size(); j++)
		{
			if (gate_vec[i]->ipt_dis[j] > SpacingDis) 
			{	
				double intpart, frapart;
				frapart = modf((double)gate_vec[i]->ipt_dis[j]/SpacingDis, &intpart);	
				int NumBufferInserted = intpart;	
			//	cout<<"NumBufferInserted: "<<NumBufferInserted<<endl;
				double SeparateDis = (double)gate_vec[i]->ipt_dis[j]/(NumBufferInserted+1);
			//	cout<<"SeparateDis: "<<SeparateDis<<endl;
				
				vector<string> BufferName;
				BufferName.resize(NumBufferInserted);
				for (int k=0; k<NumBufferInserted; k++)
				{	
					ostringstream nameinfo;
					nameinfo<<"Buffer_"<<gate_vec[i]->ipt_nodes[j]<<"_"<<k<<"_"<<gate_vec[i]->name;
					BufferName[k] = nameinfo.str();
				//	cout<<BufferName<<endl;
				
					gate *temp_gate_ptr = new gate;
					temp_gate_ptr->setGateInfo(BufferName[k],"Ibuffer",60,100);
					if (k==0) 
					{	
						temp_gate_ptr->ipt_nodes.push_back(gate_vec[i]->ipt_nodes[j]);
						temp_gate_ptr->ipt_ptrs.push_back(gate_vec[i]->ipt_ptrs[j]);
					}
					else 
					{	
						temp_gate_ptr->ipt_nodes.push_back(BufferName[k-1]);
						// the pointer to the input gate (buffer) of this current must be at the end of gate_vec since it is pushed back to the vector
						temp_gate_ptr->ipt_ptrs.push_back(gate_vec[gate_vec.size()-1]);
					}
					
					temp_gate_ptr->ipt_dis.push_back(SeparateDis);
					temp_gate_ptr->DimInit();
					gate_vec.push_back(temp_gate_ptr);
					gate_map[temp_gate_ptr->name] = temp_gate_ptr;
				}

				gate_vec[i]->ipt_dis[j] = SeparateDis;
				gate_vec[i]->ipt_nodes[j] = BufferName[NumBufferInserted-1];
				gate_vec[i]->ipt_ptrs[j] = gate_vec[gate_vec.size()-1];

			} 
		}
	}
}

/* add the pointers of fan-outs to each gate (for implementing IncrementalTimingAnalysis) */
void AddOptPtrs(vector<gate*> &gate_vec, map<string,gate*> &gate_map)
{
	for (int i=0; i<gate_vec.size(); i++)
	{
		gate* gate_ptr;
		for (int j=0; j<gate_vec[i]->ipt_ptrs.size(); j++)
		{
			gate_ptr = gate_vec[i]->ipt_ptrs[j];
			gate_ptr->opt_ptrs.push_back(gate_vec[i]);
		}
	}

	// test output
/*	for (int i=0; i<gate_vec.size(); i++)
	{
		cout<<gate_vec[i]->name<<" has output gates: ";
		gate* gate_ptr;
		for (int j=0; j<gate_vec[i]->opt_ptrs.size(); j++)
		{
			gate_ptr = gate_vec[i]->opt_ptrs[j];
			cout<<gate_ptr->name<<" ";
		}
		cout<<endl;
	}
*/
}

void TestOptimization(ofstream &outfile)
{
	double dL = 5e-9; // delta L
	double geofctr = 1.39;

	gate *cgate = gate_map0.find("n9")->second;

	if (increoption == 0)
	{	IncrGateDimArithmetic(cgate, "output", dL, MagnetLengthUB);	}
	else if (increoption == 1)
	{	IncrGateDimGeometric(cgate, "output", geofctr, MagnetLengthUB);	}
	else
	{	cout<<"Error: Undefined geofactor in TestOptimization()!!!"<<endl;	assert(false);	}

	cout<<cgate->name<<" "<<cgate->optM_L<<endl;

	vector<gate*> CriticalPathGates1;
	GateDelayCalc(gate_vec0, gate_map0);
	double MaxDelay1 = StaticTimingAnalysis(gate_vec0, gate_map0, CriticalPathGates1);
	double EnergyConsumption1 = EnergyCalc(MaxDelay1, gate_vec0);

	cout<<"TestMaxDelay: "<<MaxDelay1<<"ns"<<endl;
	cout<<"TestEnergyConsumption: "<<EnergyConsumption1<<"J"<<endl;
	cout<<"TestGates on Critical Path: ";
	
	for (int i=0; i<CriticalPathGates1.size(); i++)
	{
		cout<<CriticalPathGates1[i]->name<<" ";
	}
	cout<<endl;

	PrintGatePtrVec(gate_vec0);

	if (outfile.is_open() == 1)
	{
		outfile<<"-1 "<<MaxDelay1<<" "<<EnergyConsumption1<<endl;
	}

	if (increoption == 0)
	{	DecrGateDimArithmetic(cgate, "output", dL, MagnetLengthLB);	}
	else if (increoption == 1)
	{	DecrGateDimGeometric(cgate, "output", geofctr, MagnetLengthLB);	}
	else
	{	cout<<"Error: Undefined geofactor in TestOptimization()!!!"<<endl;	assert(false);	}
	
}

void PriOptimization(ofstream &outfile)
{
	double dL = 70e-9; // delta L
	double geofctr = 8;

	for (int i=0; i<gate_vec0.size(); i++)
	{
		if (gate_vec0[i]->gatetype == "PRI_INPUT")
		{
			if (increoption == 0)
			{			
				IncrGateDimArithmetic(gate_vec0[i], "output", dL, MagnetLengthUB);
				//	cout<<gate_vec0[i]->name<<" "<<gate_vec0[i]->optM_L<<endl;
			}
			else if (increoption == 1)
			{
				IncrGateDimGeometric(gate_vec0[i], "output", geofctr, MagnetLengthUB);
				//	cout<<gate_vec0[i]->name<<" "<<gate_vec0[i]->optM_L<<endl;
			}
			else
			{
				cout<<"Error: Undefined sizing increasing option!!!"<<endl;
				assert(false);
			}			
		}
	}

	CriticalPathGates0.clear();
	GateDelayCalc(gate_vec0, gate_map0);
	MaxDelay0 = StaticTimingAnalysis(gate_vec0, gate_map0, CriticalPathGates0);
	EnergyConsumption0 = EnergyCalc(MaxDelay0, gate_vec0);

	cout<<"PriOptimizationMaxDelay: "<<MaxDelay0<<"ns"<<endl;
	cout<<"PriOptimizationEnergyConsumption: "<<EnergyConsumption0<<"J"<<endl;
	cout<<"PriOptimizationGates on Critical Path: ";
	
	for (int i=0; i<CriticalPathGates0.size(); i++)
	{
		cout<<CriticalPathGates0[i]->name<<" ";
	}
	cout<<endl;

	if (outfile.is_open() == 1)
	{
		outfile<<"2 "<<MaxDelay0<<" "<<EnergyConsumption0<<endl;
	}


//	PrintGatePtrVec(gate_vec0);

/*	for (int i=0; i<gate_vec0.size(); i++)
	{
		if (gate_vec0[i]->gatetype == "PRI_INPUT")
		{
			if (increoption == 0)
			{
				DecrGateDimArithmetic(gate_vec0[i], "output", dL, MagnetLengthLB);
			//	cout<<gate_vec0[i]->name<<" "<<gate_vec0[i]->optM_L<<endl;
			}
			else if (increoption == 1)
			{
				DecrGateDimGeometric(gate_vec0[i], "output", dL, MagnetLengthLB);
			//	cout<<gate_vec0[i]->name<<" "<<gate_vec0[i]->optM_L<<endl;
			}
			else 
			{
				cout<<"Error: Undefined sizing decreasing option!!!"<<endl;
				assert(false);		
			} 
		}
	}
*/

}


// this is the new version, old version can be found in /Dropbox/parser for reference of the correct functionality
void Optimization(ofstream &outfile)
{
//	clock_t OPT_START, OPT_MIDDLE, OPT_END;

	int OptimizationContinueFlag = 1;
	int CURRENT_ITERATION = 0;
	int MAX_ITERATION = 40000;

//	new version
	do 
	{
	//	OPT_START = clock();

	/*	vector<double> NewDelay;
		vector<double> NewEnergyConsumption;
		vector<double> NewPower;
		vector<double> NewArea;
		NewDelay.resize(2*CriticalPathGates0.size());
		NewEnergyConsumption.resize(2*CriticalPathGates0.size());
		NewPower.resize(2*CriticalPathGates0.size());
		NewArea.resize(2*CriticalPathGates0.size());
	*/
		vector<double> TempInternalDelay;
		TempInternalDelay.resize(2*CriticalPathGates0.size());

		vector<vector<double> > TempInputInterconnectDelay;
		TempInputInterconnectDelay.resize(2*CriticalPathGates0.size());

		vector<vector<double> > TempOutputInterconnectDelay;
		TempOutputInterconnectDelay.resize(2*CriticalPathGates0.size());

/*		vector<double> SizedInternalDelay;
		SizedInternalDelay.resize(2*CriticalPathGates0.size());

		vector<vector<double> > SizedInputInterconnectDelay;
		SizedInputInterconnectDelay.resize(2*CriticalPathGates0.size());

		vector<vector<double> > SizedOutputInterconnectDelay;
		SizedOutputInterconnectDelay.resize(2*CriticalPathGates0.size());
*/
		vector<double> CriterionValue;
		CriterionValue.resize(2*CriticalPathGates0.size());

		vector< vector<gate*> > CriticalPathGatesHolder;
		CriticalPathGatesHolder.resize(2*CriticalPathGates0.size());

	//	clock_t TP1, TP2, TP3, TP4, TP5;
		
		for (int i=0; i<2*CriticalPathGates0.size(); i++)
		{
			double inpart, frapart;
			frapart = modf( (double) i/2, & inpart);

			string ioMagnetSelect;	// set the sizing quantity as input or output to a gate
			if (frapart == 0) 
			{	ioMagnetSelect = "input";	}
			else if (frapart == 0.5) 
			{	ioMagnetSelect = "output";	}
			else 
			{	assert(false);	}

			if ( (ioMagnetSelect == "input") && (std::find(gatetype1to1list.begin(), gatetype1to1list.end(), CriticalPathGates0[inpart]->gatetype) != gatetype1to1list.end()) )
			{	// if this gate is a buffer or so, only its output should be sized
				// we then set the energy/delay/CriticalPathGates value under this situation same as the unoptimized case, so it will be screened/not considered in later case
			//	NewDelay[i] = MaxDelay0;
			//	NewEnergyConsumption[i] = EnergyConsumption0;
			//	NewPower[i] = Power0;
			//	NewArea[i] = Area0;
				CriticalPathGatesHolder[i] = CriticalPathGates0;

				CriterionValue[i] = 0; //NewDelay[i];
			}
			else if ( (ioMagnetSelect == "input") && (std::find(gatetypepriinputlist.begin(), gatetypepriinputlist.end(), CriticalPathGates0[inpart]->gatetype) != gatetypepriinputlist.end()) )
			{	// if this gate is a primary input, only its output should be sized
				// we then set the energy/delay/CriticalPathGates value under this situation same as the unoptimized case, so it will be screened/not considered in later case
			//	NewDelay[i] = MaxDelay0;
			//	NewEnergyConsumption[i] = EnergyConsumption0;
			//	NewPower[i] = Power0;
			//	NewArea[i] = Area0;
				CriticalPathGatesHolder[i] = CriticalPathGates0;

				CriterionValue[i] = 0; //NewDelay[i];
			}
			else if (std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(), CriticalPathGates0[inpart]->gatetype) == gatetypeconstlist.end())
			{	// if this gate is not const (just a check), then we do the gate sizing up to see results
				
				// if IncrGateDimArithmetic/IncrGateDimGeometric function finds out the input or output magnet of this gate has already sized to its upper-bound, 
				// then its will not change the magnet size. So theoretically, the energy/delay/CriticalPathGates should not be 
				// different from unsized-up (original) case. Then by the worse case when the best delay is still equal to MaxDelay0, 
				// this one will automatically go into:
				// (NewDelay[i] == MaxDelay0) && (std::equal(CriticalPathGatesHolder[i].begin(), CriticalPathGatesHolder[i].end(), CriticalPathGates0.begin()) == 0)
				// and not considered as valid improvement. 

			//	double RelatedDelayBefore = ExtractSizedGateRelatedDelays(CriticalPathGates0[inpart], CriticalPathGates0);
			//	double RelatedPowerBefore = ExtractSizedGateRelatedPower(CriticalPathGates0[inpart]);
			//	double RelatedAreaBefore = ExtractSizedGateRelatedArea(CriticalPathGates0[inpart]);
			//	TP1 = clock();	

				vector<double> RelatedDPABefore = ExtractSizedGateRelatedDPA(CriticalPathGates0[inpart], CriticalPathGates0);

				int MaxLFlag = 0;
				if (increoption == 0)
				{	IncrGateDimArithmetic(CriticalPathGates0[inpart], ioMagnetSelect, dL, MagnetLengthUB, MaxLFlag);	}
				else if (increoption == 1)
				{	IncrGateDimGeometric(CriticalPathGates0[inpart], ioMagnetSelect, geofctr, MagnetLengthUB, MaxLFlag);	}
				else 
				{	cout<<"Error: Undefined increoption in Optimization()!!!"<<endl;	assert(false);	}
				
				if (MaxLFlag == 0) // Magnet Length Upper-bound not reached for this magnet
				{
				//	TP2 = clock();
					GateDelayExtractNCalc(CriticalPathGates0[inpart], ioMagnetSelect, TempInternalDelay[i], TempInputInterconnectDelay[i], TempOutputInterconnectDelay[i], Delay3to1TBL, Delay5to1TBL);
				//	GateDelayExtractNCalc(CriticalPathGates0[inpart], TempInternalDelay[i], TempInputInterconnectDelay[i], TempOutputInterconnectDelay[i], Delay3to1TBL, Delay5to1TBL);
				//	GateDelayExtractNCalc(CriticalPathGates0[inpart], TempInternalDelay[i], TempInputInterconnectDelay[i], TempOutputInterconnectDelay[i]);
				//	GateDelayExtractNCalc(CriticalPathGates0[inpart], TempInternalDelay[i], TempInputInterconnectDelay[i], TempOutputInterconnectDelay[i], SizedInternalDelay[i], SizedInputInterconnectDelay[i], SizedOutputInterconnectDelay[i]);
				//	GateDelayCalc(CriticalPathGates0[inpart]);	//*** new version only recalculate the delay related to the newly sized gate and go with IncrementalTA
				//	GateDelayCalc(gate_vec0, gate_map0, CriticalPathGates0[inpart]); //*** new version only recalculate the delay related to the newly sized gate
				//	TP3 = clock();

				//	GateDelayCalc(gate_vec0, gate_map0);
				//	cout<<"Delay related to sized gate after optimization: "<<CriticalPathGates0[inpart]->name<<" is "<<ExtractSizedGateRelatedDelays(CriticalPathGates0[inpart], CriticalPathGates0)<<endl;
				//	double RelatedDelayAfter = ExtractSizedGateRelatedDelays(CriticalPathGates0[inpart], CriticalPathGates0);
				//	double RelatedPowerAfter = ExtractSizedGateRelatedPower(CriticalPathGates0[inpart]);
				//	double RelatedAreaAfter = ExtractSizedGateRelatedArea(CriticalPathGates0[inpart]);
					vector<double> RelatedDPAAfter = ExtractSizedGateRelatedDPA(CriticalPathGates0[inpart], CriticalPathGates0);
				//	cout<<"Delay/power/area change related to sized gate: "<<CriticalPathGates0[inpart]->name<<" is "<<RelatedDelayAfter-RelatedDelayBefore<<" and "<<RelatedPowerAfter-RelatedPowerBefore<<" and "<<RelatedAreaAfter-RelatedAreaBefore<<endl;

				//	NewDelay[i] = StaticTimingAnalysis(gate_vec0, gate_map0, CriticalPathGatesHolder[i]);
				//	NewEnergyConsumption[i] = EnergyCalc(NewDelay[i], gate_vec0);
				//	NewPower[i] = PowerCalc(gate_vec0);
				//	NewArea[i] = AreaCalc(gate_vec0);

				//	cout<<"The actually circuit delay/power of this gate: "<<CriticalPathGates0[inpart]->name<<" is "<<NewDelay[i]<<" and "<<NewPower[i]<<" and "<<NewArea[i]<<endl;
				//	TP4 = clock();
					CriterionValue[i] = (RelatedDPAAfter[0]-RelatedDPABefore[0])/(RelatedDPAAfter[1]-RelatedDPABefore[1]); //NewDelay[i];

					if (increoption == 0)
					{	DecrGateDimArithmetic(CriticalPathGates0[inpart], ioMagnetSelect, dL, MagnetLengthLB);	}
					else if (increoption == 1)
					{	DecrGateDimGeometric(CriticalPathGates0[inpart], ioMagnetSelect, geofctr, MagnetLengthLB);	}
					else 
					{	cout<<"Error: Undefined decreoption in Optimization()!!!"<<endl;	assert(false);	}
					
					// note that every time the delay is recalculated for sized up related magnets, there has to a restorage to its original value
					// by the following function. Otherwise new delays for the circuits will be kept and not restored
					GateDelayRecover(CriticalPathGates0[inpart], TempInternalDelay[i], TempInputInterconnectDelay[i], TempOutputInterconnectDelay[i]);
				//	GateDelayCalc(CriticalPathGates0[inpart]);	//*** new version only recalculate the delay related to the newly sized gate and go with IncrementalTA
				//	GateDelayCalc(gate_vec0, gate_map0, CriticalPathGates0[inpart]); //*** new version only recalculate the delay related to the newly sized gate
				//	TP5 = clock();

				}
				else if (MaxLFlag == 1)
				{
				//	NewDelay[i] = MaxDelay0;
				//	NewEnergyConsumption[i] = EnergyConsumption0;
				//	NewPower[i] = Power0;
				//	NewArea[i] = Area0;
					CriticalPathGatesHolder[i] = CriticalPathGates0;

					CriterionValue[i] = 0; //NewDelay[i];
				}
				else
				{
					cout<<"Error: Undefined situation for Max Magnet Length Flag!!!"<<endl;
					assert(false);
				}
			}
			else 
			{
				cout<<"Error: Undefined situation in Optimization (changing gate size at the beginning section) !!!"<<endl;
				assert(false);
			}
	//		cout<<endl;
	//		cout<<"Inside attempts geo change RunTime: "<<(double)(TP2-TP1)/CLOCKS_PER_SEC<<"S"<<endl;
	//		cout<<"Inside attempts GD RunTime: "<<(double)(TP3-TP2)/CLOCKS_PER_SEC<<"S"<<endl;
	//		cout<<"Inside attempts Extract RunTime: "<<(double)(TP4-TP3)/CLOCKS_PER_SEC<<"S"<<endl;
	//		cout<<"Inside attempts Criterion and GD recover RunTime: "<<(double)(TP5-TP4)/CLOCKS_PER_SEC<<"S"<<endl;
	//		cout<<"Inside attempts a complete one RunTime: "<<(double)(TP5-TP1)/CLOCKS_PER_SEC<<"S"<<endl;
	//		cout<<endl;
		}

		// Looking for Critical Path Delay reduction
		vector<double>::iterator MinInCriterionValue = std::min_element(CriterionValue.begin(), CriterionValue.end());
		double CriterionCompareValue = 0; // MaxDelay0;

	//	OPT_MIDDLE = clock();
	//	cout<<"Attempts RunTime: "<<(double)(OPT_MIDDLE-OPT_START)/CLOCKS_PER_SEC<<"S"<<endl;

		if (*MinInCriterionValue < CriterionCompareValue) // if the delay did improves, we adapt this delay and show the dimension change
		{
			int OrderInCriterionValue = std::distance(CriterionValue.begin(), MinInCriterionValue);
			gate *ImprovedGate = CriticalPathGates0[OrderInCriterionValue/2];

			string ImprovedIO; // from the location of the improved delay in NewDelay vector, decide it is input or output
			if ( OrderInCriterionValue % 2 == 0 )
			{	
				ImprovedIO = "input"; 
			}
			else if ( OrderInCriterionValue % 2 == 1 ) 
			{
				ImprovedIO = "output"; 
			}

			if (increoption == 0)
			{	IncrGateDimArithmetic(ImprovedGate, ImprovedIO, dL, MagnetLengthUB);	}
			else if (increoption == 1)
			{	IncrGateDimGeometric(ImprovedGate, ImprovedIO, geofctr, MagnetLengthUB);	}
			else 
			{	cout<<"Error: Undefined increoption in Optimization()!!!"<<endl;	assert(false);	}
			
			cout<<"In "<<CURRENT_ITERATION<<"th Iteration, Delay Improved this time by changing: "<<ImprovedIO<<" magnet(s) of "<<ImprovedGate->name;

			if (ImprovedIO == "input")
			{
				cout<<" To: "<<ImprovedGate->iptM_L[0]<<"m"<<endl;
			}
			else if (ImprovedIO == "output")
			{
				cout<<" To: "<<ImprovedGate->optM_L<<"m"<<endl;
			}
			else 
			{
				cout<<"Error: Undefined action in Optimization function when trying to display change magnet length!!!"<<endl;
				assert(false); 
			}
		
		//	GateDelayRecover(CriticalPathGates0[OrderInCriterionValue/2], SizedInternalDelay[OrderInCriterionValue], SizedInputInterconnectDelay[OrderInCriterionValue], SizedOutputInterconnectDelay[OrderInCriterionValue]);
			GateDelayCalc(CriticalPathGates0[OrderInCriterionValue/2], Delay3to1TBL, Delay5to1TBL);
		//	GateDelayCalc(CriticalPathGates0[OrderInCriterionValue/2]);	//*** new version only recalculate the delay related to the newly sized gate and go with IncrementalTA
		//	GateDelayCalc(gate_vec0, gate_map0, CriticalPathGates0[OrderInCriterionValue/2]); //*** new version only recalculate the delay related to the newly sized gate

			CriticalPathGates0.clear(); // set new critical path
			MaxDelay0 = IncrementalTimingAnalysis(gate_vec0, gate_map0, CriticalPathGates0, CriticalPathGates0[OrderInCriterionValue/2]);
		//	MaxDelay0 = StaticTimingAnalysis(gate_vec0, gate_map0, CriticalPathGates0);

			vector<double> PowerArea = PowerAreaCalc(gate_vec0); 	// if not need power/area numebr, comment this line/function for fast speed	//	EnergyConsumption0 = EnergyCalc(MaxDelay0, gate_vec0);
	
			cout<<"In "<<CURRENT_ITERATION<<"th Iteration, Delay Improved to: "<<MaxDelay0<<"ns"<<endl;
			
			if (outfile.is_open() == 1)
			{
			//	outfile<<CURRENT_ITERATION+1<<" "<<MaxDelay0<<endl;	// if not need power/area numebr, run this line for fast speed
				outfile<<CURRENT_ITERATION+1<<" "<<MaxDelay0<<" "<<PowerArea[0]<<" "<<PowerArea[1]<<endl;
			}

		/*	cout<<"Gates on Critical Path: ";
			for (int i=0; i<CriticalPathGates0.size(); i++)
			{
				cout<<CriticalPathGates0[i]->name<<" ";
			}
			cout<<endl;
		*/
		//	PrintGatePtrVec(gate_vec0);
			OptimizationContinueFlag = 1;
			CURRENT_ITERATION++;
		}
		else if (*MinInCriterionValue >= CriterionCompareValue) // if every sizing up makes delay worse
		{
			cout<<"In "<<CURRENT_ITERATION<<"th Iteration, Delay Became Worse and No Valid Improvement Found! Try other optimization direction! OPTIMIZATION ENDS!"<<endl;	
			OptimizationContinueFlag = 0;		
		}
		else
		{
			cout<<"Error: Undefined situation inside the optimization function when judging the results! Optimization stops and program terminates!!!"<<endl;
			cout<<"MinInCriterionValue: "<<*MinInCriterionValue<<endl;
			OptimizationContinueFlag = 0;
			assert(false);
		}
		
	//	OPT_END = clock();
	//	cout<<"Adaption RunTime: "<<(double)(OPT_END-OPT_MIDDLE)/CLOCKS_PER_SEC<<"S"<<endl;

	} while ( (OptimizationContinueFlag == 1) && (CURRENT_ITERATION < MAX_ITERATION ) );

}

// Recalculation of the results with rounding
void RecalcWithRounding(ofstream &outfile)
{
	cout << endl;
	for (int i=0; i<gate_vec0.size(); i++)
	{
	//	cout << "For gate: " << gate_vec0[i]->name << ", ";
		if ( gate_vec0[i]->iptM_L.size() != 0)
		{
		//	cout << "its input magnets sizes: ";
			for (int j=0; j<gate_vec0[i]->iptM_L.size(); j++)
			{
				gate_vec0[i]->iptM_L[j] = 1e-8 * round( gate_vec0[i]->iptM_L[j] * 1e8 );
			//	cout << gate_vec0[i]->iptM_L[j] << "nm, ";			
			}
		}
		
	//	cout << "its output magnet size: ";
		gate_vec0[i]->optM_L = 1e-8 * round( gate_vec0[i]->optM_L * 1e8 );
	//	cout << gate_vec0[i]->optM_L << "nm." << endl;

	}

	GateDelayCalc(gate_vec0, gate_map0);
	MaxDelay0 = StaticTimingAnalysis(gate_vec0, gate_map0, CriticalPathGates0);
	EnergyConsumption0 = EnergyCalc(MaxDelay0, gate_vec0);
	Power0 = PowerCalc(gate_vec0);
	Area0 = AreaCalc(gate_vec0);

	cout << endl;
	cout << "Data after rounding is performed: " << endl;
	cout << "Delay after rounding: " << MaxDelay0 << "s" << endl;
	cout << "Energy after rounding: " << EnergyConsumption0 << "J" << endl;
	cout << "Power after rounding: " << Power0 << "W" << endl;
//	cout << "Area after rounding: " << Area0 << "m^2" << endl;
}

/* Operations after optimiztion*/
void PostOperation()
{
	OptimalLengthStatis(gate_vec0);
	FreePtrs(gate_vec0);
}


