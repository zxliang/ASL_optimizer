#include "gate.h"

/* Extract the delays/power/area related to the sized gate on the critical path */ // unfinished
vector<double> ExtractSizedGateRelatedDPA(gate *sgate_ptr, vector<gate*> &CriticalPathGates)
{
	// The way of power in this function is calculated following swoption = 0

	vector<double> SizedGateRelatedDPA = {0,0,0};
	double Ic_total = 0;
	double M_L_total = 0;

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
				SizedGateRelatedDPA[0] += ogate_ptr->ext_ipt_delay[position];
			}
		}

		Ic_total += Vdd/(rho_F*sgate_ptr->optM_T/(sgate_ptr->optM_W*sgate_ptr->optM_L) + rho_N*midC_T/(sgate_ptr->optM_W*sgate_ptr->optM_L));

		M_L_total += M_L_total + sgate_ptr->optM_L;
	}
	else if ( std::find(gatetype1to1list.begin(), gatetype1to1list.end(), sgate_ptr->gatetype) != gatetype1to1list.end() )
	{
		// for the delay change of the fan-in gates of this sized gate on the critical path
		SizedGateRelatedDPA[0] += sgate_ptr->ext_ipt_delay[0];

		// for the delay change of the fan-out gates of this sized gate on the critical path
		gate *ogate_ptr;
		for (int i=0; i<sgate_ptr->opt_ptrs.size(); i++)
		{
			if ( std::find(CriticalPathGates.begin(), CriticalPathGates.end(), sgate_ptr->opt_ptrs[i]) != CriticalPathGates.end() )
			{
				ogate_ptr = sgate_ptr->opt_ptrs[i];
				int position = std::find(ogate_ptr->ipt_ptrs.begin(), ogate_ptr->ipt_ptrs.end(), sgate_ptr) - ogate_ptr->ipt_ptrs.begin();
				SizedGateRelatedDPA[0] += ogate_ptr->ext_ipt_delay[position];
			}
		}

		Ic_total += Vdd/(rho_F*sgate_ptr->optM_T/(0.5*sgate_ptr->optM_W*sgate_ptr->optM_L) + rho_N*midC_T/(0.5*sgate_ptr->optM_W*sgate_ptr->optM_L));

		M_L_total += M_L_total + sgate_ptr->optM_L;		
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
				SizedGateRelatedDPA[0] += sgate_ptr->int_ipt_delay[position];
				SizedGateRelatedDPA[0] += sgate_ptr->ext_ipt_delay[position];
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
				SizedGateRelatedDPA[0] += ogate_ptr->ext_ipt_delay[position];
			}
		}

		int NumInputMagnets = 2*sgate_ptr->ipt_ptrs.size() - 1;

		Ic_total += Vdd/(rho_F*sgate_ptr->optM_T/(0.5*sgate_ptr->optM_W*sgate_ptr->optM_L) + rho_N*midC_T/(0.5*sgate_ptr->optM_W*sgate_ptr->optM_L));		
		Ic_total += NumInputMagnets*( Vdd/(rho_F*sgate_ptr->iptM_T[0]/(0.5*sgate_ptr->iptM_W[0]*sgate_ptr->iptM_L[0]) + rho_N*midC_T/(0.5*sgate_ptr->iptM_W[0]*sgate_ptr->iptM_L[0])) );

		M_L_total = M_L_total + sgate_ptr->optM_L;
		M_L_total = M_L_total + NumInputMagnets*sgate_ptr->iptM_L[0];

	}
	else
	{	
		cout<<"Error: Undefined gatetype in function ExtractSizedGateRelatedDPA!!!"<<endl;
		assert(false);
	}

	double Power = Vdd*Ic_total;
	SizedGateRelatedDPA[1] = Power;
	double Area = ioptM_W*M_L_total;
	SizedGateRelatedDPA[2] = Area;
	return SizedGateRelatedDPA;
}

double ExtractSizedGateRelatedPower(gate *sgate_ptr)
{
	double Power = 0;

	if (swoption == 0)
	{
		double Ic_total = 0;
		
		if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),sgate_ptr->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),sgate_ptr->gatetype) != gatetype5to1list.end()) ) 
		{
			Ic_total = Ic_total + Vdd/(rho_F*sgate_ptr->optM_T/(0.5*sgate_ptr->optM_W*sgate_ptr->optM_L) + rho_N*midC_T/(0.5*sgate_ptr->optM_W*sgate_ptr->optM_L));
			
			int NumInputMagnets = 2*sgate_ptr->ipt_ptrs.size() - 1;
			Ic_total = Ic_total + NumInputMagnets*( Vdd/(rho_F*sgate_ptr->iptM_T[0]/(0.5*sgate_ptr->iptM_W[0]*sgate_ptr->iptM_L[0]) + rho_N*midC_T/(0.5*sgate_ptr->iptM_W[0]*sgate_ptr->iptM_L[0])) );
		}
		else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),sgate_ptr->gatetype) != gatetype1to1list.end() )
		{
			Ic_total = Ic_total + Vdd/(rho_F*sgate_ptr->optM_T/(0.5*sgate_ptr->optM_W*sgate_ptr->optM_L) + rho_N*midC_T/(0.5*sgate_ptr->optM_W*sgate_ptr->optM_L));
		}
		else if ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),sgate_ptr->gatetype) != gatetypepriinputlist.end() )
		{
			Ic_total = Ic_total + Vdd/(rho_F*sgate_ptr->optM_T/(sgate_ptr->optM_W*sgate_ptr->optM_L) + rho_N*midC_T/(sgate_ptr->optM_W*sgate_ptr->optM_L));
		}
		else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),sgate_ptr->gatetype) != gatetypeconstlist.end() )
		{
			// do nothing about Ic_total 
		}
		else 
		{
			cout<<"Error: Gatetype in PowerCalc NOT defined!!!"<<endl;
			assert(false);
		}
	
		Power = Vdd*Ic_total;

	}
	else if (swoption == 1)
	{
		Power = 1e5;
		cout<<"Warning: ExtractSizedGateRelatedPower in under this switching option is not properly defined, set to 1e10W!!!"<<endl;
	}
	else
	{
		cout<<"Error: Switching option (swoption) defined not as 0 or 1. Please check in gate.h!!!"<<endl;
		assert(false);
	}

	return Power;
}



double ExtractSizedGateRelatedArea(gate *sgate_ptr)
{
	double M_L_total = 0;
	
	if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),sgate_ptr->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),sgate_ptr->gatetype) != gatetype5to1list.end()) ) 
	{
		M_L_total = M_L_total + sgate_ptr->optM_L;

		int NumInputMagnets = 2*sgate_ptr->ipt_ptrs.size() - 1;
		M_L_total = M_L_total + NumInputMagnets*sgate_ptr->iptM_L[0];
	}
	else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),sgate_ptr->gatetype) != gatetype1to1list.end() )
	{
		M_L_total = M_L_total + sgate_ptr->optM_L;
	}
	else if ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),sgate_ptr->gatetype) != gatetypepriinputlist.end() )
	{
		M_L_total = M_L_total + sgate_ptr->optM_L;
	}
	else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),sgate_ptr->gatetype) != gatetypeconstlist.end() )
	{
		// do nothing about Ic_total 
	}
	else 
	{
		cout<<"Error: Gatetype in AreaCalc NOT defined!!!"<<endl;
		assert(false);
	}
		
	double Area = ioptM_W*M_L_total;
	return Area;
}


double EnergyCalc(double MaxDelay, vector<gate*> &gate_vec)
{
	double EnergyConsumption = 0;


	if (swoption == 3)
	{
		double Ic_total = 0;
				
		for (int i=0; i<gate_vec.size(); i++)
		{
			if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),gate_vec[i]->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),gate_vec[i]->gatetype) != gatetype5to1list.end()) ) 
			{
				// because both the input/output magnets serve as an input magnet when considering total power/energy
				double output_R_C = rho_C*ioptC_T/(gate_vec[i]->optM_W*(input_contact_factor * gate_vec[i]->optM_L + input_contact_adjustment));
				double output_R_F = rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*(input_magnet_factor * gate_vec[i]->optM_L + input_magnet_adjustment));
				double output_R_N = rho_N*midC_T/(gate_vec[i]->optM_W*(input_nonmagnet_factor * gate_vec[i]->optM_L + input_nonmagnet_adjustment));
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
				
				int NumInputMagnets = 2*gate_vec[i]->ipt_ptrs.size() - 1;

				double input_R_C = rho_C*ioptC_T/(gate_vec[i]->iptM_W[0]*(input_contact_factor * gate_vec[i]->iptM_L[0] + input_contact_adjustment));
				double input_R_F = rho_F*gate_vec[i]->iptM_T[0]/(gate_vec[i]->iptM_W[0]*(input_magnet_factor * gate_vec[i]->iptM_L[0] + input_magnet_adjustment));
				double input_R_N = rho_N*midC_T/(gate_vec[i]->iptM_W[0]*(input_nonmagnet_factor * gate_vec[i]->iptM_L[0] + input_nonmagnet_adjustment));
				double input_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + NumInputMagnets*( Vdd/( input_R_C + input_R_F + input_R_N + input_R_G ) );

			}
			else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),gate_vec[i]->gatetype) != gatetype1to1list.end() )
			{
				double output_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);				
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
			}
			else if ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),gate_vec[i]->gatetype) != gatetypepriinputlist.end() )
			{
				double output_R_C = rho_C*ioptC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
			}
			else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),gate_vec[i]->gatetype) != gatetypeconstlist.end() )
			{
				// do nothing about Ic_total 
			}
			else 
			{
				cout<<"Error: Gatetype in PowerCalc NOT defined!!!"<<endl;
				assert(false);
			}
		}
		
		EnergyConsumption = Vdd*Ic_total*MaxDelay;

	}
	else if (swoption == 2)
	{
		double Ic_total = 0;
				
		for (int i=0; i<gate_vec.size(); i++)
		{
			if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),gate_vec[i]->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),gate_vec[i]->gatetype) != gatetype5to1list.end()) ) 
			{
				double output_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
				
				int NumInputMagnets = 2*gate_vec[i]->ipt_ptrs.size() - 1;

				double input_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]);
				double input_R_F = rho_F*gate_vec[i]->iptM_T[0]/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]);
				double input_R_N = rho_N*midC_T/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]);
				double input_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + NumInputMagnets*( Vdd/( input_R_C + input_R_F + input_R_N + input_R_G ) );

			}
			else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),gate_vec[i]->gatetype) != gatetype1to1list.end() )
			{
				double output_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);				
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
			}
			else if ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),gate_vec[i]->gatetype) != gatetypepriinputlist.end() )
			{
				double output_R_C = rho_C*ioptC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
			}
			else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),gate_vec[i]->gatetype) != gatetypeconstlist.end() )
			{
				// do nothing about Ic_total 
			}
			else 
			{
				cout<<"Error: Gatetype in PowerCalc NOT defined!!!"<<endl;
				assert(false);
			}
		}
		
		EnergyConsumption = Vdd*Ic_total*MaxDelay;

	}
	else if (swoption == 0)
	{
		double Ic_total = 0;
		
		for (int i=0; i<gate_vec.size(); i++)
		{
			if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),gate_vec[i]->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),gate_vec[i]->gatetype) != gatetype5to1list.end()) ) 
			{
				Ic_total = Ic_total + Vdd/(rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L) + rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L));
				
				int NumInputMagnets = 2*gate_vec[i]->ipt_ptrs.size() - 1;
				Ic_total = Ic_total + NumInputMagnets*( Vdd/(rho_F*gate_vec[i]->iptM_T[0]/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]) + rho_N*midC_T/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0])) );

			}
			else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),gate_vec[i]->gatetype) != gatetype1to1list.end() )
			{
				Ic_total = Ic_total + Vdd/(rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L) + rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L));
			}
			else if ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),gate_vec[i]->gatetype) != gatetypepriinputlist.end() )
			{
				Ic_total = Ic_total + Vdd/(rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L) + rho_N*midC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L));
			}
			else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),gate_vec[i]->gatetype) != gatetypeconstlist.end() )
			{
				// do nothing about Ic_total 
			}
			else 
			{
				cout<<"Error: Gatetype in EnergyCalc NOT defined!!!"<<endl;
				assert(false);
			}
		}
		
		EnergyConsumption = Vdd*Ic_total*MaxDelay;

	}
	else if (swoption == 1)
	{
		EnergyConsumption = 1e10;
		cout<<"Warning: EnergyCalc in under this switching option is not properly defined, set to 1e10J!!!"<<endl;
	}
	else
	{
		cout<<"Error: Switching option (swoption) defined not as 0 or 1. Please check in gate.h!!!"<<endl;
		assert(false);
	}

	return EnergyConsumption;
}

double PowerCalc(vector<gate*> &gate_vec)
{
	double Power = 0;


	if (swoption == 3)
	{
		double Ic_total = 0;
				
		for (int i=0; i<gate_vec.size(); i++)
		{
			if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),gate_vec[i]->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),gate_vec[i]->gatetype) != gatetype5to1list.end()) ) 
			{
				double output_R_C = rho_C*ioptC_T/(gate_vec[i]->optM_W*(input_contact_factor * gate_vec[i]->optM_L + input_contact_adjustment));
				double output_R_F = rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*(input_magnet_factor * gate_vec[i]->optM_L + input_magnet_adjustment));
				double output_R_N = rho_N*midC_T/(gate_vec[i]->optM_W*(input_nonmagnet_factor * gate_vec[i]->optM_L + input_nonmagnet_adjustment));
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
				
				int NumInputMagnets = 2*gate_vec[i]->ipt_ptrs.size() - 1;

				double input_R_C = rho_C*ioptC_T/(gate_vec[i]->iptM_W[0]*(input_contact_factor * gate_vec[i]->iptM_L[0] + input_contact_adjustment));
				double input_R_F = rho_F*gate_vec[i]->iptM_T[0]/(gate_vec[i]->iptM_W[0]*(input_magnet_factor * gate_vec[i]->iptM_L[0] + input_magnet_adjustment));
				double input_R_N = rho_N*midC_T/(gate_vec[i]->iptM_W[0]*(input_nonmagnet_factor * gate_vec[i]->iptM_L[0] + input_nonmagnet_adjustment));
				double input_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + NumInputMagnets*( Vdd/( input_R_C + input_R_F + input_R_N + input_R_G ) );

			}
			else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),gate_vec[i]->gatetype) != gatetype1to1list.end() )
			{
				double output_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);				
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
			}
			else if ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),gate_vec[i]->gatetype) != gatetypepriinputlist.end() )
			{
				double output_R_C = rho_C*ioptC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
			}
			else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),gate_vec[i]->gatetype) != gatetypeconstlist.end() )
			{
				// do nothing about Ic_total 
			}
			else 
			{
				cout<<"Error: Gatetype in PowerCalc NOT defined!!!"<<endl;
				assert(false);
			}
		}
		
		Power = Vdd*Ic_total;

	}
	else if (swoption == 2)
	{
		double Ic_total = 0;
				
		for (int i=0; i<gate_vec.size(); i++)
		{
			if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),gate_vec[i]->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),gate_vec[i]->gatetype) != gatetype5to1list.end()) ) 
			{
				double output_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
				
				int NumInputMagnets = 2*gate_vec[i]->ipt_ptrs.size() - 1;

				double input_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]);
				double input_R_F = rho_F*gate_vec[i]->iptM_T[0]/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]);
				double input_R_N = rho_N*midC_T/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]);
				double input_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + NumInputMagnets*( Vdd/( input_R_C + input_R_F + input_R_N + input_R_G ) );

			}
			else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),gate_vec[i]->gatetype) != gatetype1to1list.end() )
			{
				double output_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);				
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
			}
			else if ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),gate_vec[i]->gatetype) != gatetypepriinputlist.end() )
			{
				double output_R_C = rho_C*ioptC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total = Ic_total + Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
			}
			else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),gate_vec[i]->gatetype) != gatetypeconstlist.end() )
			{
				// do nothing about Ic_total 
			}
			else 
			{
				cout<<"Error: Gatetype in PowerCalc NOT defined!!!"<<endl;
				assert(false);
			}
		}
		
		Power = Vdd*Ic_total;

	}
	else if (swoption == 0)
	{
		double Ic_total = 0;
		
		for (int i=0; i<gate_vec.size(); i++)
		{
			if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),gate_vec[i]->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),gate_vec[i]->gatetype) != gatetype5to1list.end()) ) 
			{
				Ic_total = Ic_total + Vdd/(rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L) + rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L));
				
				int NumInputMagnets = 2*gate_vec[i]->ipt_ptrs.size() - 1;
				Ic_total = Ic_total + NumInputMagnets*( Vdd/(rho_F*gate_vec[i]->iptM_T[0]/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]) + rho_N*midC_T/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0])) );

			}
			else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),gate_vec[i]->gatetype) != gatetype1to1list.end() )
			{
				Ic_total = Ic_total + Vdd/(rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L) + rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L));
			}
			else if ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),gate_vec[i]->gatetype) != gatetypepriinputlist.end() )
			{
				Ic_total = Ic_total + Vdd/(rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L) + rho_N*midC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L));
			}
			else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),gate_vec[i]->gatetype) != gatetypeconstlist.end() )
			{
				// do nothing about Ic_total 
			}
			else 
			{
				cout<<"Error: Gatetype in PowerCalc NOT defined!!!"<<endl;
				assert(false);
			}
		}
		
		Power = Vdd*Ic_total;

	}
	else if (swoption == 1)
	{
		Power = 1e10;
		cout<<"Warning: PowerCalc in under this switching option is not properly defined, set to 1e10W!!!"<<endl;
		assert(false);
	}
	else
	{
		cout<<"Error: Switching option (swoption) defined not as 0 or 1. Please check in gate.h!!!"<<endl;
		assert(false);
	}

	return Power;
}

double AreaCalc(vector<gate*> &gate_vec)
{
	double M_L_total = 0;
	
	for (int i=0; i<gate_vec.size(); i++)
	{
		if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),gate_vec[i]->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),gate_vec[i]->gatetype) != gatetype5to1list.end()) ) 
		{
			M_L_total = M_L_total + gate_vec[i]->optM_L;
			
			int NumInputMagnets = 2*gate_vec[i]->ipt_ptrs.size() - 1;
			M_L_total = M_L_total + NumInputMagnets*gate_vec[i]->iptM_L[0];
		}
		else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),gate_vec[i]->gatetype) != gatetype1to1list.end() )
		{
			M_L_total = M_L_total + gate_vec[i]->optM_L;
		}
		else if ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),gate_vec[i]->gatetype) != gatetypepriinputlist.end() )
		{
			M_L_total = M_L_total + gate_vec[i]->optM_L;
		}
		else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),gate_vec[i]->gatetype) != gatetypeconstlist.end() )
		{
			// do nothing about Ic_total 
		}
		else 
		{
			cout<<"Error: Gatetype in AreaCalc NOT defined!!!"<<endl;
			assert(false);
		}
	}
		
	double Area = ioptM_W*M_L_total;
	return Area;
}

vector<double> PowerAreaCalc(vector<gate*> &gate_vec)
{
	double M_L_total = 0;
	double Ic_total = 0;

	if (swoption == 3)
	{
		for (int i=0; i<gate_vec.size(); i++)
		{
			if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),gate_vec[i]->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),gate_vec[i]->gatetype) != gatetype5to1list.end()) ) 
			{
				M_L_total += gate_vec[i]->optM_L;

				double output_R_C = rho_C*ioptC_T/(gate_vec[i]->optM_W*(input_contact_factor * gate_vec[i]->optM_L + input_contact_adjustment));
				double output_R_F = rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*(input_magnet_factor * gate_vec[i]->optM_L + input_magnet_adjustment));
				double output_R_N = rho_N*midC_T/(gate_vec[i]->optM_W*(input_nonmagnet_factor * gate_vec[i]->optM_L + input_nonmagnet_adjustment));
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total += Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
			
				int NumInputMagnets = 2*gate_vec[i]->ipt_ptrs.size() - 1;

				M_L_total += NumInputMagnets*gate_vec[i]->iptM_L[0];

				double input_R_C = rho_C*ioptC_T/(gate_vec[i]->iptM_W[0]*(input_contact_factor * gate_vec[i]->iptM_L[0] + input_contact_adjustment));
				double input_R_F = rho_F*gate_vec[i]->iptM_T[0]/(gate_vec[i]->iptM_W[0]*(input_magnet_factor * gate_vec[i]->iptM_L[0] + input_magnet_adjustment));
				double input_R_N = rho_N*midC_T/(gate_vec[i]->iptM_W[0]*(input_nonmagnet_factor * gate_vec[i]->iptM_L[0] + input_nonmagnet_adjustment));
				double input_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total += NumInputMagnets*( Vdd/( input_R_C + input_R_F + input_R_N + input_R_G ) );

			}
			else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),gate_vec[i]->gatetype) != gatetype1to1list.end() )
			{
				M_L_total += gate_vec[i]->optM_L;

				double output_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);				
				Ic_total += Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );

			}
			else if ( gate_vec[i]->gatetype == "PRI_INPUT" )
			{
				M_L_total += gate_vec[i]->optM_L;

				double output_R_C = rho_C*ioptC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total += Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
				
			}
			else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),gate_vec[i]->gatetype) != gatetypeconstlist.end() )
			{
				// do nothing about Ic_total 
			}
			else 
			{
				cout<<"Error: Gatetype in PowerAreaCalc NOT defined!!!"<<endl;
				assert(false);
			}
		}	
	}
	else if (swoption == 2)
	{
		for (int i=0; i<gate_vec.size(); i++)
		{
			if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),gate_vec[i]->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),gate_vec[i]->gatetype) != gatetype5to1list.end()) ) 
			{
				M_L_total += gate_vec[i]->optM_L;

				double output_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total += Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
			
				int NumInputMagnets = 2*gate_vec[i]->ipt_ptrs.size() - 1;

				M_L_total += NumInputMagnets*gate_vec[i]->iptM_L[0];

				double input_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]);
				double input_R_F = rho_F*gate_vec[i]->iptM_T[0]/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]);
				double input_R_N = rho_N*midC_T/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]);
				double input_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total += NumInputMagnets*( Vdd/( input_R_C + input_R_F + input_R_N + input_R_G ) );

			}
			else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),gate_vec[i]->gatetype) != gatetype1to1list.end() )
			{
				M_L_total += gate_vec[i]->optM_L;

				double output_R_C = rho_C*ioptC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);				
				Ic_total += Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );

			}
			else if ( gate_vec[i]->gatetype == "PRI_INPUT" )
			{
				M_L_total += gate_vec[i]->optM_L;

				double output_R_C = rho_C*ioptC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_F = rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_N = rho_N*midC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L);
				double output_R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_total += Vdd/( output_R_C + output_R_F + output_R_N + output_R_G );
				
			}
			else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),gate_vec[i]->gatetype) != gatetypeconstlist.end() )
			{
				// do nothing about Ic_total 
			}
			else 
			{
				cout<<"Error: Gatetype in PowerAreaCalc NOT defined!!!"<<endl;
				assert(false);
			}
		}		
	}
	else if (swoption == 1)
	{
		cout<<"Error: Undefined swoption == 1 in PowerAreaCalc!!!"<<endl;
		assert(false);
	}
	else if (swoption == 0)
	{
		for (int i=0; i<gate_vec.size(); i++)
		{
			if ( (std::find(gatetype3to1list.begin(),gatetype3to1list.end(),gate_vec[i]->gatetype) != gatetype3to1list.end()) || (std::find(gatetype5to1list.begin(),gatetype5to1list.end(),gate_vec[i]->gatetype) != gatetype5to1list.end()) ) 
			{
				M_L_total += gate_vec[i]->optM_L;
				Ic_total += Vdd/(rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L) + rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L));
			
				int NumInputMagnets = 2*gate_vec[i]->ipt_ptrs.size() - 1;

				M_L_total += NumInputMagnets*gate_vec[i]->iptM_L[0];
				Ic_total += NumInputMagnets*( Vdd/(rho_F*gate_vec[i]->iptM_T[0]/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0]) + rho_N*midC_T/(0.5*gate_vec[i]->iptM_W[0]*gate_vec[i]->iptM_L[0])) );
			}
			else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),gate_vec[i]->gatetype) != gatetype1to1list.end() )
			{
				M_L_total += gate_vec[i]->optM_L;
				Ic_total += Vdd/(rho_F*gate_vec[i]->optM_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L) + rho_N*midC_T/(0.5*gate_vec[i]->optM_W*gate_vec[i]->optM_L));
			}
			else if ( gate_vec[i]->gatetype == "PRI_INPUT" )
			{
				M_L_total += gate_vec[i]->optM_L;
				Ic_total += Vdd/(rho_F*gate_vec[i]->optM_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L) + rho_N*midC_T/(gate_vec[i]->optM_W*gate_vec[i]->optM_L));
			}
			else if ( std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(),gate_vec[i]->gatetype) != gatetypeconstlist.end() )
			{
				// do nothing about Ic_total 
			}
			else 
			{
				cout<<"Error: Gatetype in PowerAreaCalc NOT defined!!!"<<endl;
				assert(false);
			}
		}
	}

	double Area = ioptM_W*M_L_total;
	double Power = Vdd*Ic_total;
	vector<double> PowerArea {Power, Area};

	return PowerArea;
}




