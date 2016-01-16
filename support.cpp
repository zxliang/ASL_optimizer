# include "gate.h"

/* Count number of words in a string */
int WordCount(string s)
{
	int word_count(0);
	istringstream ss(s);
	string word;
	while (ss>>word) ++word_count;
	return word_count;
}

void PrintGateInternalSizeNDelay(vector<gate*> &gate_vec)
{
	cout<<endl;
	for (int i=0; i<gate_vec.size(); i++)
	{
		cout<<gate_vec[i]->name<<":(";
		if (gate_vec[i]->ipt_nodes.size() > 0) cout<<gate_vec[i]->iptM_L[0];
		cout<<","<<gate_vec[i]->optM_L<<") with delay:";
		if (gate_vec[i]->int_ipt_delay.size() > 0) cout<<gate_vec[i]->int_ipt_delay[0];
		cout<<endl;
	}			
	cout<<endl;
}

void PrintGatePtrVec(vector<gate*> &gate_vec)
{
	gate *gate_ptr;
	for (int i=0; i<gate_vec.size(); i++)
	{
		gate_ptr = gate_vec[i];
		cout<<gate_ptr->name<<" "<<gate_ptr->gatetype;
		cout<<" ("<<gate_ptr->optn_x<<","<<gate_ptr->optn_y<<")";
		
		if (gate_ptr->int_ipt_delay.size() != 0)
		{
			cout<<"~"<<gate_ptr->int_ipt_delay[0]<<"ns ";
		}
		
		for (int j=0; j<gate_ptr->ipt_nodes.size(); j++)
		{
		//	cout<<" "<<gate_ptr->ipt_nodes[j]; // same display as below
			cout<<" "<<(gate_ptr->ipt_ptrs[j])->name;
			cout<<"{"<<gate_ptr->ipt_dis[j]<<"}";
			cout<<"~"<<gate_ptr->ext_ipt_delay[j]<<"ns ";
		}

		cout<<endl;
	}
}

void PrintGateReachTime(vector<gate*> &gate_vec)
{
	cout<<endl;
	for (int i=0; i<gate_vec.size(); i++)
	{
		cout<<gate_vec[i]->name<<" with RTready: "<<gate_vec[i]->RTready<<" with ReachTime: "<<gate_vec[i]->ReachTime<<endl;
	}
	cout<<endl;
}

/* get a brief distribution of the final optimal legnth statistics */
void OptimalLengthStatis(vector<gate*> &gate_vec)
{
	cout<<endl;
	cout<<"The number of magnets under each length: "<<endl;

	int MagnetCount = 0;
	int Interval = 5;
	int LBvalue = 25;
	int UBvalue = 110;
	int nInterval = (UBvalue-LBvalue)/Interval;
	vector<int> LengthStatis (nInterval,0);

	gate* gate_ptr;
	for (int i=0; i<gate_vec.size(); i++)
	{
		gate_ptr = gate_vec[i];
	//	cout<<gate_vec[i]->gatetype<<endl;
		if (std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(), gate_ptr->gatetype) != gatetypeconstlist.end())
		{
			// do nothing
		//	cout<<gate_ptr->name<<" "<<gate_ptr->gatetype<<endl;
		}
		else if ( (std::find(gatetypepriinputlist.begin(), gatetypepriinputlist.end(), gate_ptr->gatetype) != gatetypepriinputlist.end()) || (std::find(gatetype1to1list.begin(), gatetype1to1list.end(), gate_ptr->gatetype) != gatetype1to1list.end()) )
		{
			MagnetCount += 1;
		//	cout<<gate_ptr->name<<" "<<gate_ptr->gatetype<<" "<<gate_ptr->iptM_L.size()<<endl;
			
			// for output magnet
			for (int n=0; n<nInterval; n++)
			{
				int LB = LBvalue+n*Interval;
				int UB = LBvalue+(n+1)*Interval;
				if ( (gate_ptr->optM_L >= LB*1e-9) && (gate_ptr->optM_L < UB*1e-9) )
				{
					LengthStatis[n] += 1;
				}
			}

						
		}
		else if ( (std::find(gatetype3to1list.begin(), gatetype3to1list.end(), gate_ptr->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(), gatetype5to1list.end(), gate_ptr->gatetype) != gatetype5to1list.end()) )
		{
			int nMagnet = 2*gate_ptr->iptM_L.size()-1;
			MagnetCount += nMagnet + 1;
		//	cout<<gate_ptr->name<<" "<<gate_ptr->gatetype<<" "<<gate_ptr->iptM_L.size()<<endl;
		// using an underlying definition that iptM_L[0] represents the length for all the input magnet length
			
			// for input magnet
			for (int n=0; n<nInterval; n++)
			{
				int LB = LBvalue+n*Interval;
				int UB = LBvalue+(n+1)*Interval;
				if ( (gate_ptr->iptM_L[0] >= LB*1e-9) && (gate_ptr->iptM_L[0] < UB*1e-9) )
				{
				//	if (gate_ptr->name == "n12") {cout<<"n12 ipt "<<LB<<" 3ge "<<UB<<endl;}
					LengthStatis[n] += nMagnet;
				}
			}


			// for output magnet
			for (int n=0; n<nInterval; n++)
			{
				int LB = LBvalue+n*Interval;
				int UB = LBvalue+(n+1)*Interval;
				if ( (gate_ptr->optM_L >= LB*1e-9) && (gate_ptr->optM_L < UB*1e-9) )
				{
				//	if (gate_ptr->name == "n12") {cout<<"n12 opt "<<LB<<" 1ge "<<UB<<endl;}
					LengthStatis[n] += 1;
				}
			}

		}
		else
		{
			cout<<"Error: undefined gatetype in OptimalLengthStatis()!!!"<<endl;
			assert(false);
		}
	}

	for (int j=0; j<nInterval; j++)
	{
		int LB = LBvalue+j*Interval;
		int UB = LBvalue+(j+1)*Interval;
		cout<<"For "<<LB<<"<=L<"<<UB<<": "<<LengthStatis[j]<<endl;
	}
	
	cout<<"Total magnet count: "<<MagnetCount<<endl;

}


/* free all the pointer and memories after the program finishes */
void FreePtrs(vector<gate*> &gate_vec)
{
	gate *gate_ptr;
	for (int i=0; i<gate_vec.size(); i++)
	{	gate_ptr = gate_vec[i];
		delete gate_ptr;
		gate_ptr = nullptr;
	}
}


 

