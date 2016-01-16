# include "gate.h"

void GCstamp(MatrixXd &G, int ipt_n, int opt_n, double C_L, double C_W, double C_T, Matrix2d &GC, Matrix2d &G0C)	// add C (for Contact) stamp to G, NOTE THAT argument order different
{
	double LC = C_T/lambda_C;
	double RC = rho_C*lambda_C/(C_W*C_L);
	
//	Matrix2d GC;
	GC(0,0) = 1/(RC*LC);
	GC(0,1) = 0;
	GC(1,0) = 0;
	GC(1,1) = 1/(RC*LC)*LC*1/sinh(LC);

//	Matrix2d G0C;
	G0C(0,0) = 0;
	G0C(0,1) = 0;
	G0C(1,0) = 0;
	G0C(1,1) = 1/RC*(1/tanh(LC) - 1/sinh(LC));

	if ((ipt_n != -1) && (opt_n != -1)){
	
		G(2*ipt_n,2*ipt_n) = G(2*ipt_n,2*ipt_n) + GC(0,0) + G0C(0,0);
		G(2*ipt_n+1,2*ipt_n) = G(2*ipt_n+1,2*ipt_n) + GC(1,0) + G0C(1,0);
		G(2*ipt_n,2*ipt_n+1) = G(2*ipt_n,2*ipt_n+1) + GC(0,1) + G0C(0,1);
		G(2*ipt_n+1,2*ipt_n+1) = G(2*ipt_n+1,2*ipt_n+1) + GC(1,1) + G0C(1,1);

		G(2*opt_n,2*opt_n) = G(2*opt_n,2*opt_n) + GC(0,0) + G0C(0,0);
		G(2*opt_n+1,2*opt_n) = G(2*opt_n+1,2*opt_n) + GC(1,0) + G0C(1,0);
		G(2*opt_n,2*opt_n+1) = G(2*opt_n,2*opt_n+1) + GC(0,1) + G0C(0,1);
		G(2*opt_n+1,2*opt_n+1) = G(2*opt_n+1,2*opt_n+1) + GC(1,1) + G0C(1,1);

		G(2*ipt_n,2*opt_n) = G(2*ipt_n,2*opt_n) - GC(0,0);
		G(2*ipt_n+1,2*opt_n) = G(2*ipt_n+1,2*opt_n) - GC(1,0);
		G(2*ipt_n,2*opt_n+1) = G(2*ipt_n,2*opt_n+1) - GC(0,1);
		G(2*ipt_n+1,2*opt_n+1) = G(2*ipt_n+1,2*opt_n+1) - GC(1,1);

		G(2*opt_n,2*ipt_n) = G(2*opt_n,2*ipt_n) - GC(0,0);
		G(2*opt_n+1,2*ipt_n) = G(2*opt_n+1,2*ipt_n) - GC(1,0);
		G(2*opt_n,2*ipt_n+1) = G(2*opt_n,2*ipt_n+1) - GC(0,1);
		G(2*opt_n+1,2*ipt_n+1) = G(2*opt_n+1,2*ipt_n+1) - GC(1,1);
	}
	else if ((ipt_n == -1) && (opt_n == -1)) {
		
		cout<<"Invalid Node Detected for Nstamp!"<<endl;
	}
	else if ((ipt_n == -1) && (opt_n != -1)) {
		
		G(2*opt_n,2*opt_n) = G(2*opt_n,2*opt_n) + GC(0,0) + G0C(0,0);
		G(2*opt_n+1,2*opt_n) = G(2*opt_n+1,2*opt_n) + GC(1,0) + G0C(1,0);
		G(2*opt_n,2*opt_n+1) = G(2*opt_n,2*opt_n+1) + GC(0,1) + G0C(0,1);
		G(2*opt_n+1,2*opt_n+1) = G(2*opt_n+1,2*opt_n+1) + GC(1,1) + G0C(1,1);
	}
	else if ((ipt_n != -1) && (opt_n == -1)) {
		
		G(2*ipt_n,2*ipt_n) = G(2*ipt_n,2*ipt_n) + GC(0,0) + G0C(0,0);
		G(2*ipt_n+1,2*ipt_n) = G(2*ipt_n+1,2*ipt_n) + GC(1,0) + G0C(1,0);
		G(2*ipt_n,2*ipt_n+1) = G(2*ipt_n,2*ipt_n+1) + GC(0,1) + G0C(0,1);
		G(2*ipt_n+1,2*ipt_n+1) = G(2*ipt_n+1,2*ipt_n+1) + GC(1,1) + G0C(1,1);
	}

}

void GFstamp(MatrixXd &G, int ipt_n, int opt_n, double M_L, double M_W, double M_T, Matrix2d &GF, Matrix2d &G0F)	// add F stamp to G
{
	double LF = M_T/lambda_F;
	double RF = rho_F*lambda_F/(M_W*M_L);
	
//	Matrix2d GF;
	GF(0,0) = 1/(RF*LF);
	GF(0,1) = 1/(RF*LF)*PF;
	GF(1,0) = 1/(RF*LF)*PF;
	GF(1,1) = 1/(RF*LF)*PF*PF + (1-PF*PF)/RF*(1/sinh(LF));

//	Matrix2d G0F;
	G0F(0,0) = 0;
	G0F(0,1) = 0;
	G0F(1,0) = 0;
	G0F(1,1) = (1-PF*PF)/RF*(1/tanh(LF) - 1/sinh(LF));

	if ((ipt_n != -1) && (opt_n != -1)) {

		G(2*ipt_n,2*ipt_n) = G(2*ipt_n,2*ipt_n) + GF(0,0) + G0F(0,0);
		G(2*ipt_n+1,2*ipt_n) = G(2*ipt_n+1,2*ipt_n) + GF(1,0) + G0F(1,0);
		G(2*ipt_n,2*ipt_n+1) = G(2*ipt_n,2*ipt_n+1) + GF(0,1) + G0F(0,1);
		G(2*ipt_n+1,2*ipt_n+1) = G(2*ipt_n+1,2*ipt_n+1) + GF(1,1) + G0F(1,1);

		G(2*opt_n,2*opt_n) = G(2*opt_n,2*opt_n) + GF(0,0) + G0F(0,0);
		G(2*opt_n+1,2*opt_n) = G(2*opt_n+1,2*opt_n) + GF(1,0) + G0F(1,0);
		G(2*opt_n,2*opt_n+1) = G(2*opt_n,2*opt_n+1) + GF(0,1) + G0F(0,1);
		G(2*opt_n+1,2*opt_n+1) = G(2*opt_n+1,2*opt_n+1) + GF(1,1) + G0F(1,1);
	
		G(2*ipt_n,2*opt_n) = G(2*ipt_n,2*opt_n) - GF(0,0);
		G(2*ipt_n+1,2*opt_n) = G(2*ipt_n+1,2*opt_n) - GF(1,0);
		G(2*ipt_n,2*opt_n+1) = G(2*ipt_n,2*opt_n+1) - GF(0,1);
		G(2*ipt_n+1,2*opt_n+1) = G(2*ipt_n+1,2*opt_n+1) - GF(1,1);
	
		G(2*opt_n,2*ipt_n) = G(2*opt_n,2*ipt_n) - GF(0,0);
		G(2*opt_n+1,2*ipt_n) = G(2*opt_n+1,2*ipt_n) - GF(1,0);
		G(2*opt_n,2*ipt_n+1) = G(2*opt_n,2*ipt_n+1) - GF(0,1);
		G(2*opt_n+1,2*ipt_n+1) = G(2*opt_n+1,2*ipt_n+1) - GF(1,1);
	}
	else if ((ipt_n == -1) && (opt_n == -1)) {
		
		cout<<"Invalid Node Detected for Fstamp!"<<endl;
	}
	else if ((ipt_n == -1) && (opt_n != -1)) {
		
		G(2*opt_n,2*opt_n) = G(2*opt_n,2*opt_n) + GF(0,0) + G0F(0,0);
		G(2*opt_n+1,2*opt_n) = G(2*opt_n+1,2*opt_n) + GF(1,0) + G0F(1,0);
		G(2*opt_n,2*opt_n+1) = G(2*opt_n,2*opt_n+1) + GF(0,1) + G0F(0,1);
		G(2*opt_n+1,2*opt_n+1) = G(2*opt_n+1,2*opt_n+1) + GF(1,1) + G0F(1,1);
	}
	else if ((ipt_n != -1) && (opt_n == -1)) {
		
		G(2*ipt_n,2*ipt_n) = G(2*ipt_n,2*ipt_n) + GF(0,0) + G0F(0,0);
		G(2*ipt_n+1,2*ipt_n) = G(2*ipt_n+1,2*ipt_n) + GF(1,0) + G0F(1,0);
		G(2*ipt_n,2*ipt_n+1) = G(2*ipt_n,2*ipt_n+1) + GF(0,1) + G0F(0,1);
		G(2*ipt_n+1,2*ipt_n+1) = G(2*ipt_n+1,2*ipt_n+1) + GF(1,1) + G0F(1,1);
	}

}

void GNstamp(MatrixXd &G, int ipt_n, int opt_n, double N_W, double N_T, double N_L, Matrix2d &GN, Matrix2d &G0N)	// add N stamp to G, NOTE THAT argument order different
{
	double LN = N_L/lambda_N;
	double RN = rho_N*lambda_N/(N_W*N_T);
	
//	Matrix2d GN;
	GN(0,0) = 1/(RN*LN);
	GN(0,1) = 0;
	GN(1,0) = 0;
	GN(1,1) = 1/(RN*LN)*LN*1/sinh(LN);

//	Matrix2d G0N;
	G0N(0,0) = 0;
	G0N(0,1) = 0;
	G0N(1,0) = 0;
	G0N(1,1) = 1/RN*(1/tanh(LN) - 1/sinh(LN));

	if ((ipt_n != -1) && (opt_n != -1)){
	
		G(2*ipt_n,2*ipt_n) = G(2*ipt_n,2*ipt_n) + GN(0,0) + G0N(0,0);
		G(2*ipt_n+1,2*ipt_n) = G(2*ipt_n+1,2*ipt_n) + GN(1,0) + G0N(1,0);
		G(2*ipt_n,2*ipt_n+1) = G(2*ipt_n,2*ipt_n+1) + GN(0,1) + G0N(0,1);
		G(2*ipt_n+1,2*ipt_n+1) = G(2*ipt_n+1,2*ipt_n+1) + GN(1,1) + G0N(1,1);

		G(2*opt_n,2*opt_n) = G(2*opt_n,2*opt_n) + GN(0,0) + G0N(0,0);
		G(2*opt_n+1,2*opt_n) = G(2*opt_n+1,2*opt_n) + GN(1,0) + G0N(1,0);
		G(2*opt_n,2*opt_n+1) = G(2*opt_n,2*opt_n+1) + GN(0,1) + G0N(0,1);
		G(2*opt_n+1,2*opt_n+1) = G(2*opt_n+1,2*opt_n+1) + GN(1,1) + G0N(1,1);

		G(2*ipt_n,2*opt_n) = G(2*ipt_n,2*opt_n) - GN(0,0);
		G(2*ipt_n+1,2*opt_n) = G(2*ipt_n+1,2*opt_n) - GN(1,0);
		G(2*ipt_n,2*opt_n+1) = G(2*ipt_n,2*opt_n+1) - GN(0,1);
		G(2*ipt_n+1,2*opt_n+1) = G(2*ipt_n+1,2*opt_n+1) - GN(1,1);

		G(2*opt_n,2*ipt_n) = G(2*opt_n,2*ipt_n) - GN(0,0);
		G(2*opt_n+1,2*ipt_n) = G(2*opt_n+1,2*ipt_n) - GN(1,0);
		G(2*opt_n,2*ipt_n+1) = G(2*opt_n,2*ipt_n+1) - GN(0,1);
		G(2*opt_n+1,2*ipt_n+1) = G(2*opt_n+1,2*ipt_n+1) - GN(1,1);
	}
	else if ((ipt_n == -1) && (opt_n == -1)) {
		
		cout<<"Invalid Node Detected for Nstamp!"<<endl;
	}
	else if ((ipt_n == -1) && (opt_n != -1)) {
		
		G(2*opt_n,2*opt_n) = G(2*opt_n,2*opt_n) + GN(0,0) + G0N(0,0);
		G(2*opt_n+1,2*opt_n) = G(2*opt_n+1,2*opt_n) + GN(1,0) + G0N(1,0);
		G(2*opt_n,2*opt_n+1) = G(2*opt_n,2*opt_n+1) + GN(0,1) + G0N(0,1);
		G(2*opt_n+1,2*opt_n+1) = G(2*opt_n+1,2*opt_n+1) + GN(1,1) + G0N(1,1);
	}
	else if ((ipt_n != -1) && (opt_n == -1)) {
		
		G(2*ipt_n,2*ipt_n) = G(2*ipt_n,2*ipt_n) + GN(0,0) + G0N(0,0);
		G(2*ipt_n+1,2*ipt_n) = G(2*ipt_n+1,2*ipt_n) + GN(1,0) + G0N(1,0);
		G(2*ipt_n,2*ipt_n+1) = G(2*ipt_n,2*ipt_n+1) + GN(0,1) + G0N(0,1);
		G(2*ipt_n+1,2*ipt_n+1) = G(2*ipt_n+1,2*ipt_n+1) + GN(1,1) + G0N(1,1);
	}

}

void GGstamp(MatrixXd &G, int ipt_n, int opt_n, double G_W, double G_T, double G_L, Matrix2d &GG, Matrix2d &G0G)	// add G (for Ground) stamp to G, NOTE THAT argument order different
{
	double LG = G_L/lambda_G;
	double RG = rho_G*lambda_G/(G_W*G_T);
	
//	Matrix2d GG;
	GG(0,0) = 1/(RG*LG);
	GG(0,1) = 0;
	GG(1,0) = 0;
	GG(1,1) = 1/(RG*LG)*LG*1/sinh(LG);

//	Matrix2d G0G;
	G0G(0,0) = 0;
	G0G(0,1) = 0;
	G0G(1,0) = 0;
	G0G(1,1) = 1/RG*(1/tanh(LG) - 1/sinh(LG));

	if ((ipt_n != -1) && (opt_n != -1)){
	
		G(2*ipt_n,2*ipt_n) = G(2*ipt_n,2*ipt_n) + GG(0,0) + G0G(0,0);
		G(2*ipt_n+1,2*ipt_n) = G(2*ipt_n+1,2*ipt_n) + GG(1,0) + G0G(1,0);
		G(2*ipt_n,2*ipt_n+1) = G(2*ipt_n,2*ipt_n+1) + GG(0,1) + G0G(0,1);
		G(2*ipt_n+1,2*ipt_n+1) = G(2*ipt_n+1,2*ipt_n+1) + GG(1,1) + G0G(1,1);

		G(2*opt_n,2*opt_n) = G(2*opt_n,2*opt_n) + GG(0,0) + G0G(0,0);
		G(2*opt_n+1,2*opt_n) = G(2*opt_n+1,2*opt_n) + GG(1,0) + G0G(1,0);
		G(2*opt_n,2*opt_n+1) = G(2*opt_n,2*opt_n+1) + GG(0,1) + G0G(0,1);
		G(2*opt_n+1,2*opt_n+1) = G(2*opt_n+1,2*opt_n+1) + GG(1,1) + G0G(1,1);

		G(2*ipt_n,2*opt_n) = G(2*ipt_n,2*opt_n) - GG(0,0);
		G(2*ipt_n+1,2*opt_n) = G(2*ipt_n+1,2*opt_n) - GG(1,0);
		G(2*ipt_n,2*opt_n+1) = G(2*ipt_n,2*opt_n+1) - GG(0,1);
		G(2*ipt_n+1,2*opt_n+1) = G(2*ipt_n+1,2*opt_n+1) - GG(1,1);

		G(2*opt_n,2*ipt_n) = G(2*opt_n,2*ipt_n) - GG(0,0);
		G(2*opt_n+1,2*ipt_n) = G(2*opt_n+1,2*ipt_n) - GG(1,0);
		G(2*opt_n,2*ipt_n+1) = G(2*opt_n,2*ipt_n+1) - GG(0,1);
		G(2*opt_n+1,2*ipt_n+1) = G(2*opt_n+1,2*ipt_n+1) - GG(1,1);
	}
	else if ((ipt_n == -1) && (opt_n == -1)) {
		
		cout<<"Invalid Node Detected for Nstamp!"<<endl;
	}
	else if ((ipt_n == -1) && (opt_n != -1)) {
		
		G(2*opt_n,2*opt_n) = G(2*opt_n,2*opt_n) + GG(0,0) + G0G(0,0);
		G(2*opt_n+1,2*opt_n) = G(2*opt_n+1,2*opt_n) + GG(1,0) + G0G(1,0);
		G(2*opt_n,2*opt_n+1) = G(2*opt_n,2*opt_n+1) + GG(0,1) + G0G(0,1);
		G(2*opt_n+1,2*opt_n+1) = G(2*opt_n+1,2*opt_n+1) + GG(1,1) + G0G(1,1);
	}
	else if ((ipt_n != -1) && (opt_n == -1)) {
		
		G(2*ipt_n,2*ipt_n) = G(2*ipt_n,2*ipt_n) + GG(0,0) + G0G(0,0);
		G(2*ipt_n+1,2*ipt_n) = G(2*ipt_n+1,2*ipt_n) + GG(1,0) + G0G(1,0);
		G(2*ipt_n,2*ipt_n+1) = G(2*ipt_n,2*ipt_n+1) + GG(0,1) + G0G(0,1);
		G(2*ipt_n+1,2*ipt_n+1) = G(2*ipt_n+1,2*ipt_n+1) + GG(1,1) + G0G(1,1);
	}

}
void Vstamp(MatrixXd &G, int V_n, int n)	// n means how many nodes before this voltage node is inserted, so it should increase 
{
	G(2*V_n,2*n) = G(2*V_n,2*n) + 1;
	G(2*V_n+1,2*n+1) = G(2*V_n+1,2*n+1) + 1;

	G(2*n,2*V_n) = G(2*n,2*V_n) + 1;
	G(2*n+1,2*V_n+1) = G(2*n+1,2*V_n+1) + 1;

}

MatrixXd ReshapeV(VectorXd &V)	// reshape the result vector into a 2xN (N according to the type of the gate)
{
//	cout<<V.size()<<endl;
	int n_col = 0.5*V.size();
	
	MatrixXd R_V = MatrixXd::Zero(2,n_col);
	for (int i=0; i<n_col; i++) {
		R_V(0,i) = V(2*i);
		R_V(1,i) = V(2*i+1);
	}
	
	return R_V;
}



