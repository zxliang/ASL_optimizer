#include "gate.h"

// First three functions are arithmetic factor sizing vesion
// version without an MaxLengthFlag 
void IncrGateDimArithmetic(gate *gate_ptr, string ioselect, double DeltaL, double MagnetLengthUpperBound)
{
	if (std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(), gate_ptr->gatetype) != gatetypeconstlist.end())
	{
		cout<<"Error: These constant type should not be considered in magnet length increase/decrease in IncrGateDimArithmetic function!!!"<<endl;
		assert(false);
	}
	else if ( (ioselect == "input") && (std::find(gatetype1to1list.begin(), gatetype1to1list.end(), gate_ptr->gatetype) == gatetype1to1list.end()) && (std::find(gatetypepriinputlist.begin(), gatetypepriinputlist.end(), gate_ptr->gatetype) == gatetypepriinputlist.end()) )
	{
		for (int i=0; i<gate_ptr->iptM_L.size(); i++)
		{
			if (gate_ptr->iptM_L[i] + DeltaL > MagnetLengthUpperBound)
			{	
				cout<<"Magnet Length Upper-Bound reached for input of gate: "<<gate_ptr->name<<": "<<gate_ptr->iptM_L[i]<<", will not be increased any more!"<<endl;
			}
			else 
			{
				gate_ptr->iptM_L[i] = gate_ptr->iptM_L[i] + DeltaL;
			}
		}
	}
	else if (ioselect == "output")
	{
		if (gate_ptr->optM_L + DeltaL > MagnetLengthUpperBound)
		{	
			cout<<"Magnet Length Upper-Bound reached for output of gate: "<<gate_ptr->name<<": "<<gate_ptr->optM_L<<", will not be increased any more!"<<endl;
		}
		else 
		{
			gate_ptr->optM_L = gate_ptr->optM_L + DeltaL;
		}	
	}
	else 
	{	
		cout<<"Error: Not defined situation for gate "<<gate_ptr->name<<" in IncrGateDimArithmetic function!!!"<<endl;
		assert(false);
	}
}

// version with an MaxLengthFlag
void IncrGateDimArithmetic(gate *gate_ptr, string ioselect, double DeltaL, double MagnetLengthUpperBound, int &MaxLengthFlag)
{
	MaxLengthFlag = 0;

	if (std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(), gate_ptr->gatetype) != gatetypeconstlist.end())
	{
		cout<<"Error: These constant/primary_input type should not be considered in magnet length increase/decrease in IncrGateDimArithmetic function!!!"<<endl;
		assert(false);
	}
	else if ( (ioselect == "input") && (std::find(gatetype1to1list.begin(), gatetype1to1list.end(), gate_ptr->gatetype) == gatetype1to1list.end()) && (std::find(gatetypepriinputlist.begin(), gatetypepriinputlist.end(), gate_ptr->gatetype) == gatetypepriinputlist.end()) )
	{
		for (int i=0; i<gate_ptr->iptM_L.size(); i++)
		{
			if (gate_ptr->iptM_L[i] + DeltaL > MagnetLengthUpperBound)
			{	
				cout<<"Magnet Length Upper-Bound reached for input of gate: "<<gate_ptr->name<<": "<<gate_ptr->iptM_L[i]<<", will not be increased any more!"<<endl;
				MaxLengthFlag = 1;
			}
			else 
			{
				gate_ptr->iptM_L[i] = gate_ptr->iptM_L[i] + DeltaL;
			}
		}
	}
	else if (ioselect == "output")
	{
		if (gate_ptr->optM_L + DeltaL > MagnetLengthUpperBound)
		{	
			cout<<"Magnet Length Upper-Bound reached for output of gate: "<<gate_ptr->name<<": "<<gate_ptr->optM_L<<", will not be increased any more!"<<endl;
			MaxLengthFlag = 1;
		}
		else 
		{
			gate_ptr->optM_L = gate_ptr->optM_L + DeltaL;
		}	
	}
	else 
	{	
		cout<<"Error: Not defined situation for gate "<<gate_ptr->name<<" in IncrGateDimArithmetic function!!!"<<endl;
		assert(false);
	}
}

void DecrGateDimArithmetic(gate *gate_ptr, string ioselect, double DeltaL, double MagnetLengthLowerBound)
{
	if (std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(), gate_ptr->gatetype) != gatetypeconstlist.end())
	{
		cout<<"Error: These constant/primary_input type should not be considered in magnet length increase/decrease in DecrGateDimArithmetic function!!!"<<endl;
		assert(false);
	}
	else if ( (ioselect == "input") && (std::find(gatetype1to1list.begin(), gatetype1to1list.end(), gate_ptr->gatetype) == gatetype1to1list.end()) && (std::find(gatetypepriinputlist.begin(), gatetypepriinputlist.end(), gate_ptr->gatetype) == gatetypepriinputlist.end()) )
	{
		for (int i=0; i<gate_ptr->iptM_L.size(); i++)
		{
			if (gate_ptr->iptM_L[i] - DeltaL < MagnetLengthLowerBound)
			{	
				cout<<"Magnet Length Lower-Bound reached for input of gate: "<<gate_ptr->name<<": "<<gate_ptr->iptM_L[i]<<", will not be decreased any more!"<<endl;
			}
			else 
			{
				gate_ptr->iptM_L[i] = gate_ptr->iptM_L[i] - DeltaL;
			}
		}
	}
	else if (ioselect == "output")
	{
		if (gate_ptr->optM_L - DeltaL < MagnetLengthLowerBound)
		{	
			cout<<"Magnet Length Lower-Bound reached for output of gate: "<<gate_ptr->name<<": "<<gate_ptr->optM_L<<", will not be decreased any more!"<<endl;
		}
		else 
		{
			gate_ptr->optM_L = gate_ptr->optM_L - DeltaL;
		}	
	}
	else 
	{	
		cout<<"Error: Not defined situation for gate "<<gate_ptr->name<<" in DecrGateDimArithmetic function!!!"<<endl;
		assert(false);
	}

}

// netxt three functions are geometric factor sizing version
// version without an MaxLengthFlag 
void IncrGateDimGeometric(gate *gate_ptr, string ioselect, double geofactor, double MagnetLengthUpperBound)
{
	if (std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(), gate_ptr->gatetype) != gatetypeconstlist.end())
	{
		cout<<"Error: These constant type should not be considered in magnet length increase/decrease in IncrGateDimGeometic function!!!"<<endl;
		assert(false);
	}
	else if ( (ioselect == "input") && (std::find(gatetype1to1list.begin(), gatetype1to1list.end(), gate_ptr->gatetype) == gatetype1to1list.end()) && (std::find(gatetypepriinputlist.begin(), gatetypepriinputlist.end(), gate_ptr->gatetype) == gatetypepriinputlist.end()) )
	{
		for (int i=0; i<gate_ptr->iptM_L.size(); i++)
		{
			if (gate_ptr->iptM_L[i] * geofactor > MagnetLengthUpperBound)
			{	
				cout<<"Magnet Length Upper-Bound reached for input of gate: "<<gate_ptr->name<<": "<<gate_ptr->iptM_L[i]<<", will not be increased any more!"<<endl;
			}
			else 
			{
				gate_ptr->iptM_L[i] = gate_ptr->iptM_L[i] * geofactor;
			}
		}
	}
	else if (ioselect == "output")
	{
		if (gate_ptr->optM_L * geofactor > MagnetLengthUpperBound)
		{	
			cout<<"Magnet Length Upper-Bound reached for output of gate: "<<gate_ptr->name<<": "<<gate_ptr->optM_L<<", will not be increased any more!"<<endl;
		}
		else 
		{
			gate_ptr->optM_L = gate_ptr->optM_L * geofactor;
		}	
	}
	else 
	{	
		cout<<"Error: Not defined situation for gate "<<gate_ptr->name<<" in IncrGateDimGeometic function!!!"<<endl;
		assert(false);
	}
}

// version with an MaxLengthFlag
void IncrGateDimGeometric(gate *gate_ptr, string ioselect, double geofactor, double MagnetLengthUpperBound, int &MaxLengthFlag)
{
	MaxLengthFlag = 0;

	if (std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(), gate_ptr->gatetype) != gatetypeconstlist.end())
	{
		cout<<"Error: These constant/primary_input type should not be considered in magnet length increase/decrease in IncrGateDimGeometric function!!!"<<endl;
		assert(false);
	}
	else if ( (ioselect == "input") && (std::find(gatetype1to1list.begin(), gatetype1to1list.end(), gate_ptr->gatetype) == gatetype1to1list.end()) && (std::find(gatetypepriinputlist.begin(), gatetypepriinputlist.end(), gate_ptr->gatetype) == gatetypepriinputlist.end()) )
	{
		for (int i=0; i<gate_ptr->iptM_L.size(); i++)
		{
			if (gate_ptr->iptM_L[i] * geofactor > MagnetLengthUpperBound)
			{	
				cout<<"Magnet Length Upper-Bound reached for input of gate: "<<gate_ptr->name<<": "<<gate_ptr->iptM_L[i]<<", will not be increased any more!"<<endl;
				MaxLengthFlag = 1;
			}
			else 
			{
				gate_ptr->iptM_L[i] = gate_ptr->iptM_L[i] * geofactor;
			}
		}
	}
	else if (ioselect == "output")
	{
		if (gate_ptr->optM_L * geofactor > MagnetLengthUpperBound)
		{	
			cout<<"Magnet Length Upper-Bound reached for output of gate: "<<gate_ptr->name<<": "<<gate_ptr->optM_L<<", will not be increased any more!"<<endl;
			MaxLengthFlag = 1;
		}
		else 
		{
			gate_ptr->optM_L = gate_ptr->optM_L * geofactor;
		}	
	}
	else 
	{	
		cout<<"Error: Not defined situation for gate "<<gate_ptr->name<<" in IncrGateDimGeometric function!!!"<<endl;
		assert(false);
	}
}

void DecrGateDimGeometric(gate *gate_ptr, string ioselect, double geofactor, double MagnetLengthLowerBound)
{
	if (std::find(gatetypeconstlist.begin(), gatetypeconstlist.end(), gate_ptr->gatetype) != gatetypeconstlist.end())
	{
		cout<<"Error: These constant/primary_input type should not be considered in magnet length increase/decrease in DecrGateDimGeometric function!!!"<<endl;
		assert(false);
	}
	else if ( (ioselect == "input") && (std::find(gatetype1to1list.begin(), gatetype1to1list.end(), gate_ptr->gatetype) == gatetype1to1list.end()) && (std::find(gatetypepriinputlist.begin(), gatetypepriinputlist.end(), gate_ptr->gatetype) == gatetypepriinputlist.end()) )
	{
		for (int i=0; i<gate_ptr->iptM_L.size(); i++)
		{
			if (gate_ptr->iptM_L[i] / geofactor < MagnetLengthLowerBound)
			{	
				cout<<"Magnet Length Lower-Bound reached for input of gate: "<<gate_ptr->name<<": "<<gate_ptr->iptM_L[i]<<", will not be decreased any more!"<<endl;
			}
			else 
			{
				gate_ptr->iptM_L[i] = gate_ptr->iptM_L[i] / geofactor;
			}
		}
	}
	else if (ioselect == "output")
	{
		if (gate_ptr->optM_L / geofactor < MagnetLengthLowerBound)
		{	
			cout<<"Magnet Length Lower-Bound reached for output of gate: "<<gate_ptr->name<<": "<<gate_ptr->optM_L<<", will not be decreased any more!"<<endl;
		}
		else 
		{
			gate_ptr->optM_L = gate_ptr->optM_L / geofactor;
		}	
	}
	else 
	{	
		cout<<"Error: Not defined situation for gate "<<gate_ptr->name<<" in DecrGateDimGeometric function!!!"<<endl;
		assert(false);
	}

}
