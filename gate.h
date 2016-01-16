#ifndef GATE_H
#define GATE_H
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <algorithm>

using namespace std;

// Eigen added 03/22/2015 for matrix computation
#include <Eigen/Dense>
#include <Eigen/LU>
using namespace Eigen;

// global constants defined (physical and material)
const double q = 1.6e-19;	// [C] electron
const double Ms = 780e+3;	// [A/m] or [emu/cc] (1 emu/cc = 1000 A/m) saturation magnetization)
const double mu_B = 9.274e-24;	// [J/T] Bohr magneton

const double PF = 0.6;	// polarization factor for magnets

const double rho_C = 7e-9;	// resistivity of contacts
const double rho_F = 170e-9;	// resistivity of ferromagnets
const double rho_N = 7e-9;	// resistivity of non-magnet channels
const double rho_G = 7e-9;	// resistivity of grounds

const double lambda_C = 400e-9;	// spin-diffusion length of contacts
const double lambda_F = 5e-9;	// spin-diffusion length of ferromagnets
const double lambda_N = 400e-9;	// spin-diffusion length of non-magnet channels
const double lambda_G = 400e-9;	// spin-diffusion length of grounds

// global constants (fixed dimensional parameters)
const double ioptC_W = 10e-9;	// globally same width of input/output contacts
const double ioptC_L = 30e-9;	// globally same length of input/output contacts (initialization?)
const double ioptC_T = 5e-9;	// globally same thickness of input/output contacts

const double ioptM_W = 10e-9;	// globally same width of input/ouput magnets
const double ioptM_L = 30e-9;	// globally same length of input/output magnets for convenient initialization
const double ioptM_T = 5e-9;	// globally same thichness of input/output magnets

const double midC_W = 20e-9;	// globally same width of middle channel connecting input/output magnets
const double midC_L = 55e-9;	// globally same length of middle channel connecting input/output magnets (only used for initialization)
const double midC_T = 30e-9;	// globally same thickness of middle channel connecting input/output magnets

const double iptG_W = 140e-9;	// globally same width of grounds at input side
const double iptG_L = 2800e-9;	// globally same length of grounds at input side for L/lambda_G consideration
const double iptG_T = 140e-9;	// globally same thichness of grounds at input side

// global constants (simulation related)
const double Vdd = 0.03;	// global voltage
const double Ic = 0.01;	// current source 
const double fsw = 4.8;	// additional switching factor

// global switching mechanism setting
const int swoption = 3;	// 1 use voltage source, 0 use current source, 2 use Icc+voltage based spin injection efficiency

// global control factor for input/output magnet portion
const double input_magnet_factor = 1;	// 0.5 for swoption == 2, 1 for swoption == 3
const double output_magnet_factor = 0.5;	// 0.5 for both swoption == 2/3
const double input_magnet_deduction = 0;	// only for swoption == 2
const double output_magnet_deduction = 0;	// only for swoption == 2
const double input_magnet_adjustment = -10e-9;	// only for swoption == 3
const double output_magnet_adjustment = 5e-9;	// only for swoption == 3

const double input_contact_factor = 1;	// used only for swoption == 3
const double output_contact_factor = 1;	// used only for swoption == 3
const double input_contact_adjustment = 0;	// only for swoption == 3
const double output_contact_adjustment = 0;	// only for swoption == 3

const double input_nonmagnet_factor = 1;	// used only for swoption == 3
const double output_nonmagnet_factor = 0;	// used only for swoption == 3
const double input_nonmagnet_adjustment = -20e-9;	// only for swoption == 3
const double output_nonmagnet_adjustment = 10e-9;	// only for swoption == 3

// global length incremental setting
const int increoption = 1; // 1 use geometric factor, 0 use arithmetic  factor

// global geometry optimization paramters
const double dL = 5e-9; // delta L
const double geofctr = 1.37; // the geometric factor to multiply in every iteration

// gate type definition for global use
const vector<string> gatetype3to1list{"na2a0","na2b0","na2c0","na2d0","no2a0","no2b0","no2c0","no2d0","no2e0"};
const vector<string> gatetype5to1list{"na3a0","na3b0","na3c0","na3d0","na3e0"};
// "Ibuffer" for inserted buffer
const vector<string> gatetype1to1list{"bufa0","bufb0","bufc0","bufd0","bufe0","nota0","notb0","notc0","notd0","note0","Ibuffer"};
const vector<string> gatetypeconstlist{"ZERO","ONE"};
const vector<string> gatetypepriinputlist{"PRI_INPUT"};

// Upper-Bound for the sizing of the length of magnets
const double MagnetLengthUB = 110e-9;
const double MagnetLengthLB = 29e-9;

class gate{
public:
	string name;
	string gatetype;
	int width, height;
	int optn_x, optn_y;
		
	vector<string> ipt_nodes;
	vector<gate*> ipt_ptrs;
	vector<gate*> opt_ptrs;

	vector<int> ipt_loc_x, ipt_loc_y;	// location of the input gates to this gate
	vector<int> ipt_dis;		// distance of input gates to this gate

	double optM_L, optM_W, optM_T; // Length/Width/Thickness of the output magnet
	vector<double> iptM_L, iptM_W, iptM_T; // Length/Width/Thickness of the input magnets
	vector<double> iptC_L, iptC_W, iptC_T; // Length/Width/Thickness of the input channels

	vector<double> InputReachTime; // record the reach time from each input, set the largest one as reachtime

	// delay is calculated after the definition of magnet/channel dimensions
	vector<double> int_ipt_delay, ext_ipt_delay; // the delay for each input inside (internal int) the gate and from earlier input (external, ext) gates 
	
	// note that both gate and primary input have these two values
	// a double set as the arrive time to this gate (initial as the int_ipt_delay[0] in SetGateMap()
	double ReachTime;
	// a indicator for whether this gates' reach time is correct set or not
	int RTready;

	gate(string n="unknown", string t="unknown", int w=60, int h=100):name(n), gatetype(t), width(w), height(h){}
	
	void setGateInfo(string n, string t, int w, int h){
		name=n;
		gatetype=t;
		width=w;
		height=h;
	}

	void DimInit()	// initialize dimensions of output magnet, input magnets and channels; this initialization function must be run AFTER input gates identified
	{
		if ( (ipt_nodes.size() == 0) && (gatetype != "PRI_INPUT") )
		{	
			cout<<"ipt_nodes NOT SET! INVALID DIMENSION INITIALIZATION!"<<endl;
			assert(false);
		}
		else if ( (ipt_nodes.size() == 0) && (gatetype == "PRI_INPUT") )
		{	
			// output magnet
			optM_W = ioptM_W;
			optM_L = ioptM_L;
			optM_T = ioptM_T;
		}
		else if (ipt_nodes.size() != 0 )
		{
			// output magnet
			optM_W = ioptM_W;
			optM_L = ioptM_L;
			optM_T = ioptM_T;

			// input magnets
			iptM_W.resize(ipt_nodes.size());	fill(iptM_W.begin(),iptM_W.end(),ioptM_W);
			iptM_L.resize(ipt_nodes.size());	fill(iptM_L.begin(),iptM_L.end(),ioptM_L);
			iptM_T.resize(ipt_nodes.size());	fill(iptM_T.begin(),iptM_T.end(),ioptM_T);

			// channels
			iptC_W.resize(ipt_nodes.size());	fill(iptC_W.begin(),iptC_W.end(),midC_W);
			iptC_L.resize(ipt_nodes.size());	fill(iptC_L.begin(),iptC_L.end(),midC_L);
			iptC_T.resize(ipt_nodes.size());	fill(iptC_T.begin(),iptC_T.end(),midC_T);
		}
	}

};

/* mylib.cpp */
void Read_Blif_N_Pl(ifstream &bliff, ifstream &plf);
void FindLocation(ifstream &plf, map<string, gate*> &gate_map);
void DistanceCalc(vector<gate*> &gate_vec, map<string, gate*> &gate_map);

void Initialization(ofstream &outfile);
void BufferInsertion(vector<gate*> &gate_vec, map<string,gate*> &gate_map, double SpacingDis);
void AddOptPtrs(vector<gate*> &gate_vec, map<string,gate*> &gate_map);

void TestOptimization(ofstream &outfile);
void PriOptimization(ofstream &outfile);
void Optimization(ofstream &outfile);

void RecalcWithRounding(ofstream &outfile);
void PostOperation();

/* sizing.cpp */
void IncrGateDimArithmetic(gate *gate_ptr, string ioselect, double DeltaL, double MagnetLengthUpperBound);
void IncrGateDimArithmetic(gate *gate_ptr, string ioselect, double DeltaL, double MagnetLengthUpperBound, int &MaxLengthFlag);
void DecrGateDimArithmetic(gate *gate_ptr, string ioselect, double DeltaL, double MagnetLengthLowerBound);

void IncrGateDimGeometric(gate *gate_ptr, string ioselect, double geofactor, double MagnetLengthUpperBound);
void IncrGateDimGeometric(gate *gate_ptr, string ioselect, double geofactor, double MagnetLengthUpperBound, int &MaxLengthFlag);
void DecrGateDimGeometric(gate *gate_ptr, string ioselect, double geofactor, double MagnetLengthLowerBound);

/* STA.cpp */
double IncrementalTimingAnalysis(vector<gate*> &gate_vec, map<string, gate*> &gate_map, vector<gate*> &CriticalPathGates, gate *sgate_ptr);
void ReachTimeReadyReset(gate *gate_ptr);

double StaticTimingAnalysis(vector<gate*> &gate_vec, map<string, gate*> &gate_map, vector<gate*> &CriticalPathGates);
double FindReachTime(gate *gate_ptr);
void FindCriticalPathGates(gate *gate_ptr, vector<gate*> &CriticalPathGates);

/* EPA.cpp (Energy/Power/Area related functions */
vector<double> ExtractSizedGateRelatedDPA(gate *sgate_ptr, vector<gate*> &CriticalPathGates);

double ExtractSizedGateRelatedPower(gate *sgate_ptr);
double ExtractSizedGateRelatedArea(gate *sgate_ptr);

double EnergyCalc(double MaxDelay, vector<gate*> &gate_vec);
double PowerCalc(vector<gate*> &gate_vec);
double AreaCalc(vector<gate*> &gate_vec);

vector<double> PowerAreaCalc(vector<gate*> &gate_vec);

/* delay.cpp */
void GateDelayCalc(vector<gate*> &gate_vec, map<string, gate*> &gate_map);

// BE CAREFUL: The following function can only be used when initialization is finished
// their are used since we only perform delay calculate for the gate whose dimension is changed
void GateDelayCalc(vector<gate*> &gate_vec, map<string, gate*> &gate_map, gate *sgate_ptr);
void GateDelayCalc(gate *sgate_ptr);
void GateDelayCalc(gate *sgate_ptr, map<pair<int, int>, double> &Delay3to1Table, map<pair<int, int>, double> &Delay5to1Table);

void GateDelayExtractNCalc(gate *sgate_ptr, double &InternalDelay, vector<double> &InputDelay, vector<double> &OutputDelay);
void GateDelayExtractNCalc(gate *sgate_ptr, double &InternalDelay, vector<double> &InputDelay, vector<double> &OutputDelay, map<pair<int, int>, double> &Delay3to1Table, map<pair<int, int>, double> &Delay5to1Table);
void GateDelayExtractNCalc(gate *sgate_ptr, string ioMagnetSelect, double &InternalDelay, vector<double> &InputDelay, vector<double> &OutputDelay, map<pair<int, int>, double> &Delay3to1Table, map<pair<int, int>, double> &Delay5to1Table);

void GateDelayExtractNCalc(gate *sgate_ptr, double &TempInternalDelay, vector<double> &TempInputDelay, vector<double> &TempOutputDelay, double &SizedInternalDelay, vector<double> &SizedInputDelay, vector<double> &SizedOutputDelay);
void GateDelayRecover(gate *sgate_ptr, double &TempInternalDelay, vector<double> &TempInputDelay, vector<double> &TempOutputDelay);

double ExtractSizedGateRelatedDelays(gate *sgate_ptr, vector<gate*> &CriticalPathGates);

/* delaycalc.cpp */
void InternalDelayPreCalc(double geofactor, map<pair<int, int>, double> &Delay3to1Table, map<pair<int, int>, double> &Delay5to1Table, double MagnetLengthUpperBound);
void IntDelayCalcMatch(gate *igate_ptr, map<pair<int, int>, double> &Delay3to1Table, map<pair<int, int>, double> &Delay5to1Table);

void IntDelayCalc(gate *igate_ptr);
void InterconnectDelayCalc(gate *igate_ptr, gate *ogate_ptr, int ipt_order);

/* MNA.cpp */
void GCstamp(MatrixXd &G, int ipt_n, int opt_n, double C_L, double C_W, double C_T, Matrix2d &GC, Matrix2d &G0C);
void GFstamp(MatrixXd &G, int ipt_n, int opt_n, double M_L, double M_W, double M_T, Matrix2d &GF, Matrix2d &G0F);
void GNstamp(MatrixXd &G, int ipt_n, int opt_n, double N_W, double N_T, double N_L, Matrix2d &GN, Matrix2d &G0N);
void GGstamp(MatrixXd &G, int ipt_n, int opt_n, double G_W, double G_T, double G_L, Matrix2d &GG, Matrix2d &G0G);
void Vstamp(MatrixXd &G, int V_n, int n);
MatrixXd ReshapeV(VectorXd &V);

/* support.cpp */
int WordCount(string s);
void PrintGateInternalSizeNDelay(vector<gate*> &gate_vec);
void PrintGatePtrVec(vector<gate*> &gate_vec);
void PrintGateReachTime(vector<gate*> &gate_vec);

void OptimalLengthStatis(vector<gate*> &gate_vec);
void FreePtrs(vector<gate*> &gate_vec);


#endif
