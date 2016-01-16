# include "gate.h"

void InternalDelayPreCalc(double geofactor, map<pair<int, int>, double> &Delay3to1Table, map<pair<int, int>, double> &Delay5to1Table, double MagnetLengthUpperBound)
{
	double InputMagnetLength, OutputMagnetLength;
	double t_sw3to1, t_sw5to1;
	int i, j;

	i = 0;
	while( ioptM_L*pow(geofactor, i) <= MagnetLengthUpperBound ) 
	{
		InputMagnetLength = ioptM_L*pow(geofactor, i);
		j = 0;
		while( ioptM_L*pow(geofactor, j) <= MagnetLengthUpperBound ) 
		{
			OutputMagnetLength = ioptM_L*pow(geofactor, j);
	
			pair<int, int> CurrentPair ((int)(InputMagnetLength*1e9), (int)(OutputMagnetLength*1e9));

			// delay of 3 to 1 gate
			if (swoption == 3)	// modeling 2 includes contacts on every magnet and grounds under input side 
			{	// rules to write the netlist same as in swoption 0 or 1

				int n=13;
				MatrixXd G = MatrixXd::Zero(2*(n+4),2*(n+4));
				// **** for now, we use the magnet dimensions and global parameters to define contact dimensions for simplicity of codes ****
				Matrix2d GC0,G0C0,GC1,G0C1,GC2,G0C2,GC3,G0C3;	// contacts
				GCstamp(G,0,1,input_contact_factor * InputMagnetLength + input_contact_adjustment,ioptM_W,ioptC_T,GC0,G0C0);
				GCstamp(G,3,4,input_contact_factor * InputMagnetLength + input_contact_adjustment,ioptM_W,ioptC_T,GC1,G0C1);	
				GCstamp(G,6,7,input_contact_factor * InputMagnetLength + input_contact_adjustment,ioptM_W,ioptC_T,GC2,G0C2);
				GCstamp(G,9,10,output_contact_factor * OutputMagnetLength + output_contact_adjustment,ioptM_W,ioptC_T,GC3,G0C3);

				Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3;	// magnets
				GFstamp(G,1,2,input_magnet_factor * InputMagnetLength + input_magnet_adjustment,ioptM_W,ioptM_T,GF0,G0F0);
				GFstamp(G,4,5,input_magnet_factor * InputMagnetLength + input_magnet_adjustment,ioptM_W,ioptM_T,GF1,G0F1);	
				GFstamp(G,7,8,input_magnet_factor * InputMagnetLength + input_magnet_adjustment,ioptM_W,ioptM_T,GF2,G0F2);
				GFstamp(G,10,11,output_magnet_factor * OutputMagnetLength + output_magnet_adjustment,ioptM_W,ioptM_T,GF3,G0F3);

				Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3;	// channels
				GNstamp(G,2,12,midC_W,midC_T,midC_L,GN0,G0N0);
				GNstamp(G,5,12,midC_W,midC_T,midC_L,GN1,G0N1);
				GNstamp(G,8,12,midC_W,midC_T,midC_L,GN2,G0N2);
				GNstamp(G,11,12,midC_W,midC_T,midC_L,GN3,G0N3);
				// **** for now, we use global parameters to define ground dimensions for simplicity of codes ****
				Matrix2d GN4,G0N4,GN5,G0N5,GN6,G0N6;	// grounds
				GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN4,G0N4);
				GGstamp(G,5,-1,iptG_W,iptG_T,iptG_L,GN5,G0N5);
				GGstamp(G,8,-1,iptG_W,iptG_T,iptG_L,GN6,G0N6);			

				Vstamp(G,0,13);
				Vstamp(G,3,14);
				Vstamp(G,6,15);
				Vstamp(G,9,16);

				VectorXd b(2*(n+4));
				b << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     Vdd, 0, -Vdd, 0, Vdd, 0, 0, 0;

				VectorXd x_V = G.inverse() * b;	
			//	cout<<G<<endl;
			//	cout<<x_V<<endl;
				MatrixXd r_V = ReshapeV(x_V);	
			//	cout<<r_V<<endl;	

				Vector2d I1, I7, I10;
				I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
				I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

				I7(0) = GF2(0,0)*(r_V(0,7)-r_V(0,8)) + GF2(0,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(0,0)*r_V(0,8) - G0F2(0,1)*r_V(1,8);
				I7(1) = GF2(1,0)*(r_V(0,7)-r_V(0,8)) + GF2(1,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(1,0)*r_V(0,8) - G0F2(1,1)*r_V(1,8);
			
				I10(0) = GF3(0,0)*(r_V(0,11)-r_V(0,10)) + GF3(0,1)*(r_V(1,11)-r_V(1,10)); // - G0F3(0,0)*r_V(0,10) - G0F3(0,1)*r_V(1,10);
				I10(1) = GF3(1,0)*(r_V(0,11)-r_V(0,10)) + GF3(1,1)*(r_V(1,11)-r_V(1,10)); // - G0F3(1,0)*r_V(0,10) - G0F3(1,1)*r_V(1,10);			

			//	cout<<I10(1)<<endl;
				
				double SIE = I10(1)/I1(0);
			//	cout <<SIE<<endl;

				double R_C = rho_C*ioptC_T/(ioptM_W*(input_contact_factor * InputMagnetLength + input_contact_adjustment));	// contact resistance considered
				double R_F = rho_F*ioptM_T/(ioptM_W*(input_magnet_factor * InputMagnetLength + input_magnet_adjustment));	// ferromagnet resistance considered
				double R_N = rho_N*midC_T/(ioptM_W*(input_nonmagnet_factor * InputMagnetLength + input_nonmagnet_adjustment));	// channel resistance considered
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);	// ground resistance considered
				double Ic_real = Vdd/( R_C + R_F + R_N + R_G );

				double Ns = Ms *OutputMagnetLength *ioptM_W *ioptM_T /mu_B;
				t_sw3to1 = 2*fsw*q*Ns/abs(Ic_real*SIE);
			}
			else if (swoption == 2)	// modeling 2 includes contacts on every magnet and grounds under input side 
			{	// rules to write the netlist same as in swoption 0 or 1

				int n=13;
				MatrixXd G = MatrixXd::Zero(2*(n+4),2*(n+4));
				// **** for now, we use the magnet dimensions and global parameters to define contact dimensions for simplicity of codes ****
				Matrix2d GC0,G0C0,GC1,G0C1,GC2,G0C2,GC3,G0C3;	// contacts
				GCstamp(G,0,1,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptC_T,GC0,G0C0);
				GCstamp(G,3,4,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptC_T,GC1,G0C1);	
				GCstamp(G,6,7,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptC_T,GC2,G0C2);
				GCstamp(G,9,10,output_magnet_factor * OutputMagnetLength - output_magnet_deduction,ioptM_W,ioptC_T,GC3,G0C3);

				Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3;	// magnets
				GFstamp(G,1,2,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptM_T,GF0,G0F0);
				GFstamp(G,4,5,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptM_T,GF1,G0F1);	
				GFstamp(G,7,8,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptM_T,GF2,G0F2);
				GFstamp(G,10,11,output_magnet_factor * OutputMagnetLength - output_magnet_deduction,ioptM_W,ioptM_T,GF3,G0F3);

				Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3;	// channels
				GNstamp(G,2,12,midC_W,midC_T,midC_L,GN0,G0N0);
				GNstamp(G,5,12,midC_W,midC_T,midC_L,GN1,G0N1);
				GNstamp(G,8,12,midC_W,midC_T,midC_L,GN2,G0N2);
				GNstamp(G,11,12,midC_W,midC_T,midC_L,GN3,G0N3);
				// **** for now, we use global parameters to define ground dimensions for simplicity of codes ****
				Matrix2d GN4,G0N4,GN5,G0N5,GN6,G0N6;	// grounds
				GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN4,G0N4);
				GGstamp(G,5,-1,iptG_W,iptG_T,iptG_L,GN5,G0N5);
				GGstamp(G,8,-1,iptG_W,iptG_T,iptG_L,GN6,G0N6);			

				Vstamp(G,0,13);
				Vstamp(G,3,14);
				Vstamp(G,6,15);
				Vstamp(G,9,16);

				VectorXd b(2*(n+4));
				b << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     Vdd, 0, -Vdd, 0, Vdd, 0, 0, 0;

				VectorXd x_V = G.inverse() * b;	
			//	cout<<G<<endl;
			//	cout<<x_V<<endl;
				MatrixXd r_V = ReshapeV(x_V);	
			//	cout<<r_V<<endl;	

				Vector2d I1, I7, I10;
				I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
				I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

				I7(0) = GF2(0,0)*(r_V(0,7)-r_V(0,8)) + GF2(0,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(0,0)*r_V(0,8) - G0F2(0,1)*r_V(1,8);
				I7(1) = GF2(1,0)*(r_V(0,7)-r_V(0,8)) + GF2(1,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(1,0)*r_V(0,8) - G0F2(1,1)*r_V(1,8);
			
				I10(0) = GF3(0,0)*(r_V(0,11)-r_V(0,10)) + GF3(0,1)*(r_V(1,11)-r_V(1,10)); // - G0F3(0,0)*r_V(0,10) - G0F3(0,1)*r_V(1,10);
				I10(1) = GF3(1,0)*(r_V(0,11)-r_V(0,10)) + GF3(1,1)*(r_V(1,11)-r_V(1,10)); // - G0F3(1,0)*r_V(0,10) - G0F3(1,1)*r_V(1,10);			

			//	cout<<I10(1)<<endl;
				
				double SIE = I10(1)/I1(0);
			//	cout <<SIE<<endl;

				double R_C = rho_C*ioptC_T/(0.5*ioptM_W*InputMagnetLength);	// contact resistance considered
				double R_F = rho_F*ioptM_T/(0.5*ioptM_W*InputMagnetLength);	// ferromagnet resistance considered
				double R_N = rho_N*midC_T/(0.5*ioptM_W*InputMagnetLength);	// channel resistance considered
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);	// ground resistance considered
				double Ic_real = Vdd/( R_C + R_F + R_N + R_G );

				double Ns = Ms *OutputMagnetLength *ioptM_W *ioptM_T /mu_B;
				t_sw3to1 = 2*fsw*q*Ns/abs(Ic_real*SIE);
			}
			else if (swoption == 1) 
			{	// depending on we choose to impose constant current source or voltage source, the matrices will be different
				int n=10;	// # of nodes
				MatrixXd G = MatrixXd::Zero(2*(n+3),2*(n+3));	// initialize with zeros every time, otherwise it seems it will accumulate??

				Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3;
				GFstamp(G,0,1,InputMagnetLength,ioptM_W,ioptM_T,GF0,G0F0);
				GFstamp(G,2,3,InputMagnetLength,ioptM_W,ioptM_T,GF1,G0F1);	// set the size of fixed magnet same as the first input magnet
				GFstamp(G,4,5,InputMagnetLength,ioptM_W,ioptM_T,GF2,G0F2);
				GFstamp(G,6,7,OutputMagnetLength,ioptM_W,ioptM_T,GF3,G0F3);
			
				Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3;
				GNstamp(G,1,8,midC_W,midC_T,midC_L,GN0,G0N0);
				GNstamp(G,3,8,midC_W,midC_T,midC_L,GN1,G0N1);	// set the size of input chanenl of fixed magnet same as the channel of first input magnet
				GNstamp(G,5,8,midC_W,midC_T,midC_L,GN2,G0N2);
				GNstamp(G,7,8,midC_W,midC_T,midC_L,GN3,G0N3);	// the channel length after joint set as 50e-9
			
				Matrix2d GN4,G0N4,GN5,G0N5,GN6,G0N6,GN7,G0N7;
				GNstamp(G,1,-1,midC_W,midC_T,100*lambda_N,GN4,G0N4);
				GNstamp(G,3,-1,midC_W,midC_T,100*lambda_N,GN5,G0N5);
				GNstamp(G,5,-1,midC_W,midC_T,100*lambda_N,GN6,G0N6);
				GNstamp(G,7,9,midC_W,midC_T,100*lambda_N,GN7,G0N7); 

				Vstamp(G,0,10);	// two input numbers, first is the ID/order of the node attached to the voltage source, second is how many nodes already there in the G matrix
				Vstamp(G,2,11);
				Vstamp(G,4,12);

				VectorXd b(2*(n+3));
				b << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     Vdd, 0, Vdd, 0, Vdd, 0;

				VectorXd x_V = G.inverse() * b;	
			//	cout<<G<<endl;
			//	cout<<x_V<<endl;
				MatrixXd r_V = ReshapeV(x_V);	
			//	cout<<r_V<<endl;	

			
				Vector2d I0, I4, I6;
				I0(0) = GF0(0,0)*(r_V(0,0)-r_V(0,1)) + GF0(0,1)*(r_V(1,0)-r_V(1,1)) - G0F0(0,0)*r_V(0,1) - G0F0(0,1)*r_V(1,1);
				I0(1) = GF0(1,0)*(r_V(0,0)-r_V(0,1)) + GF0(1,1)*(r_V(1,0)-r_V(1,1)) - G0F0(1,0)*r_V(0,1) - G0F0(1,1)*r_V(1,1);

				I4(0) = GF2(0,0)*(r_V(0,4)-r_V(0,5)) + GF2(0,1)*(r_V(1,4)-r_V(1,5)) - G0F2(0,0)*r_V(0,5) - G0F2(0,1)*r_V(1,5);
				I4(1) = GF2(1,0)*(r_V(0,4)-r_V(0,5)) + GF2(1,1)*(r_V(1,4)-r_V(1,5)) - G0F2(1,0)*r_V(0,5) - G0F2(1,1)*r_V(1,5);
			
				I6(0) = GF3(0,0)*(r_V(0,6)-r_V(0,7)) + GF3(0,1)*(r_V(1,6)-r_V(1,7)) - G0F3(0,0)*r_V(0,7) - G0F3(0,1)*r_V(1,7);
				I6(1) = GF3(1,0)*(r_V(0,6)-r_V(0,7)) + GF3(1,1)*(r_V(1,6)-r_V(1,7)) - G0F3(1,0)*r_V(0,7) - G0F3(1,1)*r_V(1,7);
		
			//	cout<<I6(1)<<endl;
				
				double SIE = I6(1)/I0(0);
			//	cout <<SIE<<endl;
			
				double Ns = Ms *OutputMagnetLength *ioptM_W *ioptM_T /mu_B;
				t_sw3to1 = 2*fsw*q*Ns/abs(I6(1));

			}
			else if (swoption == 0) 
			{
				int n=10;	// # of nodes
				MatrixXd G = MatrixXd::Zero(2*n,2*n);	// initialize with zeros every time, otherwise it seems it will accumulate??

				Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3;
				GFstamp(G,0,1,InputMagnetLength,ioptM_W,ioptM_T,GF0,G0F0);
				GFstamp(G,2,3,InputMagnetLength,ioptM_W,ioptM_T,GF1,G0F1);	// set the size of fixed magnet same as the first input magnet
				GFstamp(G,4,5,InputMagnetLength,ioptM_W,ioptM_T,GF2,G0F2);
				GFstamp(G,6,7,OutputMagnetLength,ioptM_W,ioptM_T,GF3,G0F3);
			
				Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3;
				GNstamp(G,1,8,midC_W,midC_T,midC_L,GN0,G0N0);
				GNstamp(G,3,8,midC_W,midC_T,midC_L,GN1,G0N1);	// set the size of input chanenl of fixed magnet same as the channel of first input magnet
				GNstamp(G,5,8,midC_W,midC_T,midC_L,GN2,G0N2);
				GNstamp(G,7,8,midC_W,midC_T,midC_L,GN3,G0N3);	// the channel length after joint set as 50e-9
			
				Matrix2d GN4,G0N4,GN5,G0N5,GN6,G0N6,GN7,G0N7;
				GNstamp(G,1,-1,midC_W,midC_T,100*lambda_N,GN4,G0N4);
				GNstamp(G,3,-1,midC_W,midC_T,100*lambda_N,GN5,G0N5);
				GNstamp(G,5,-1,midC_W,midC_T,100*lambda_N,GN6,G0N6);
				GNstamp(G,7,9,midC_W,midC_T,100*lambda_N,GN7,G0N7); 

				VectorXd C(2*n);
				C << Ic, 0, 0, 0, -Ic, 0, 0, 0, Ic, 0, 0, 0, 	// worse case delay should be considered, therefore two positve Ic and one negative Ic is used here
				       0, 0, 0, 0, 0, 0, 0, 0;

				VectorXd x_V = G.inverse() * C;	
			//	cout<<G<<endl;
			//	cout<<x_V<<endl;
				MatrixXd r_V = ReshapeV(x_V);	
			//	cout<<r_V<<endl;	

				Vector2d I0, I4, I6;
				I0(0) = GF0(0,0)*(r_V(0,0)-r_V(0,1)) + GF0(0,1)*(r_V(1,0)-r_V(1,1)) - G0F0(0,0)*r_V(0,1) - G0F0(0,1)*r_V(1,1);
				I0(1) = GF0(1,0)*(r_V(0,0)-r_V(0,1)) + GF0(1,1)*(r_V(1,0)-r_V(1,1)) - G0F0(1,0)*r_V(0,1) - G0F0(1,1)*r_V(1,1);

				I4(0) = GF2(0,0)*(r_V(0,4)-r_V(0,5)) + GF2(0,1)*(r_V(1,4)-r_V(1,5)) - G0F2(0,0)*r_V(0,5) - G0F2(0,1)*r_V(1,5);
				I4(1) = GF2(1,0)*(r_V(0,4)-r_V(0,5)) + GF2(1,1)*(r_V(1,4)-r_V(1,5)) - G0F2(1,0)*r_V(0,5) - G0F2(1,1)*r_V(1,5);
			
				I6(0) = GF3(0,0)*(r_V(0,6)-r_V(0,7)) + GF3(0,1)*(r_V(1,6)-r_V(1,7)) - G0F3(0,0)*r_V(0,7) - G0F3(0,1)*r_V(1,7);
				I6(1) = GF3(1,0)*(r_V(0,6)-r_V(0,7)) + GF3(1,1)*(r_V(1,6)-r_V(1,7)) - G0F3(1,0)*r_V(0,7) - G0F3(1,1)*r_V(1,7);
		
				double SIE = I6(1)/I0(0);
				// half of the area for the magnet and channel under the source is considered		
				double Ic_real = Vdd/(rho_F*ioptM_T/(0.5*ioptM_W*InputMagnetLength) + rho_N*midC_T/(0.5*ioptM_W*InputMagnetLength));

				double Ns = Ms *OutputMagnetLength *ioptM_W *ioptM_T /mu_B;
				t_sw3to1 = 2*fsw*q*Ns/abs(Ic_real*SIE);

			}

			// delay of 5 to 1 gate
			if (swoption == 3)
			{
				int n=19;
				MatrixXd G = MatrixXd::Zero(2*(n+6),2*(n+6));
				// **** for now, we use the magnet dimensions and global parameters to define contact dimensions for simplicity of codes ****
				Matrix2d GC0,G0C0,GC1,G0C1,GC2,G0C2,GC3,G0C3,GC4,G0C4,GC5,G0C5;	// contacts
				GCstamp(G,0,1,input_contact_factor * InputMagnetLength + input_contact_adjustment,ioptM_W,ioptC_T,GC0,G0C0);
				GCstamp(G,3,4,input_contact_factor * InputMagnetLength + input_contact_adjustment,ioptM_W,ioptC_T,GC1,G0C1);	
				GCstamp(G,6,7,input_contact_factor * InputMagnetLength + input_contact_adjustment,ioptM_W,ioptC_T,GC2,G0C2);
				GCstamp(G,9,10,input_contact_factor * InputMagnetLength + input_contact_adjustment,ioptM_W,ioptC_T,GC3,G0C3);	
				GCstamp(G,12,13,input_contact_factor * InputMagnetLength + input_contact_adjustment,ioptM_W,ioptC_T,GC4,G0C4);
				GCstamp(G,15,16,output_contact_factor * OutputMagnetLength + output_contact_adjustment,ioptM_W,ioptC_T,GC5,G0C5);

				Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3,GF4,G0F4,GF5,G0F5;	// magnets
				GFstamp(G,1,2,input_magnet_factor * InputMagnetLength + input_magnet_adjustment,ioptM_W,ioptM_T,GF0,G0F0);
				GFstamp(G,4,5,input_magnet_factor * InputMagnetLength + input_magnet_adjustment,ioptM_W,ioptM_T,GF1,G0F1);	
				GFstamp(G,7,8,input_magnet_factor * InputMagnetLength + input_magnet_adjustment,ioptM_W,ioptM_T,GF2,G0F2);
				GFstamp(G,10,11,input_magnet_factor * InputMagnetLength + input_magnet_adjustment,ioptM_W,ioptM_T,GF3,G0F3);	
				GFstamp(G,13,14,input_magnet_factor * InputMagnetLength + input_magnet_adjustment,ioptM_W,ioptM_T,GF4,G0F4);
				GFstamp(G,16,17,output_magnet_factor * OutputMagnetLength + output_magnet_adjustment,ioptM_W,ioptM_T,GF5,G0F5);

				Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3,GN4,G0N4,GN5,G0N5;	// channels
				GNstamp(G,2,18,midC_W,midC_T,midC_L,GN0,G0N0);
				GNstamp(G,5,18,midC_W,midC_T,midC_L,GN1,G0N1);
				GNstamp(G,8,18,midC_W,midC_T,midC_L,GN2,G0N2);
				GNstamp(G,11,18,midC_W,midC_T,midC_L,GN3,G0N3);
				GNstamp(G,14,18,midC_W,midC_T,midC_L,GN4,G0N4);
				GNstamp(G,17,18,midC_W,midC_T,midC_L,GN5,G0N5);
				// **** for now, we use global parameters to define ground dimensions for simplicity of codes ****
				Matrix2d GN6,G0N6,GN7,G0N7,GN8,G0N8,GN9,G0N9,GN10,G0N10;	// grounds
				GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN6,G0N6);
				GGstamp(G,5,-1,iptG_W,iptG_T,iptG_L,GN7,G0N7);
				GGstamp(G,8,-1,iptG_W,iptG_T,iptG_L,GN8,G0N8);
				GGstamp(G,11,-1,iptG_W,iptG_T,iptG_L,GN9,G0N9);
				GGstamp(G,14,-1,iptG_W,iptG_T,iptG_L,GN10,G0N10);

				Vstamp(G,0,19);
				Vstamp(G,3,20);
				Vstamp(G,6,21);
				Vstamp(G,9,22);
				Vstamp(G,12,23);
				Vstamp(G,15,24);

				VectorXd b(2*(n+6));
				b << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     Vdd, 0, -Vdd, 0, Vdd, 0, -Vdd, 0, Vdd, 0, 0, 0;

				VectorXd x_V = G.inverse() * b;	
			//	cout<<G<<endl;
			//	cout<<x_V<<endl;
				MatrixXd r_V = ReshapeV(x_V);	
			//	cout<<r_V<<endl;	

				Vector2d I1, I7, I13, I16;
				I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
				I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

				I7(0) = GF2(0,0)*(r_V(0,7)-r_V(0,8)) + GF2(0,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(0,0)*r_V(0,8) - G0F2(0,1)*r_V(1,8);
				I7(1) = GF2(1,0)*(r_V(0,7)-r_V(0,8)) + GF2(1,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(1,0)*r_V(0,8) - G0F2(1,1)*r_V(1,8);

				I13(0) = GF4(0,0)*(r_V(0,13)-r_V(0,14)) + GF4(0,1)*(r_V(1,13)-r_V(1,14)); // - G0F4(0,0)*r_V(0,13) - G0F4(0,1)*r_V(1,14);
				I13(1) = GF4(1,0)*(r_V(0,13)-r_V(0,14)) + GF4(1,1)*(r_V(1,13)-r_V(1,14)); // - G0F4(1,0)*r_V(0,13) - G0F4(1,1)*r_V(1,14);
			
				I16(0) = GF5(0,0)*(r_V(0,17)-r_V(0,16)) + GF5(0,1)*(r_V(1,17)-r_V(1,16)); // - G0F5(0,0)*r_V(0,16) - G0F5(0,1)*r_V(1,16);
				I16(1) = GF5(1,0)*(r_V(0,17)-r_V(0,16)) + GF5(1,1)*(r_V(1,17)-r_V(1,16)); // - G0F5(1,0)*r_V(0,16) - G0F5(1,1)*r_V(1,16);			

			//	cout<<I10(1)<<endl;
				
				double SIE = I16(1)/I1(0);
			//	cout <<SIE<<endl;

				double R_C = rho_C*ioptC_T/(ioptM_W*(input_contact_factor * InputMagnetLength + input_contact_adjustment));	// contact resistance considered
				double R_F = rho_F*ioptM_T/(ioptM_W*(input_magnet_factor * InputMagnetLength + input_magnet_adjustment));	// ferromagnet resistance considered
				double R_N = rho_N*midC_T/(ioptM_W*(input_nonmagnet_factor * InputMagnetLength + input_nonmagnet_adjustment));	// channel resistance considered
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);	// ground resistance considered
				double Ic_real = Vdd/( R_C + R_F + R_N + R_G );

				double Ns = Ms*OutputMagnetLength*ioptM_W*ioptM_T/mu_B;
				t_sw5to1 = 2*fsw*q*Ns/abs(Ic_real*SIE);

			}
			else if (swoption == 2)
			{
				int n=19;
				MatrixXd G = MatrixXd::Zero(2*(n+6),2*(n+6));
				// **** for now, we use the magnet dimensions and global parameters to define contact dimensions for simplicity of codes ****
				Matrix2d GC0,G0C0,GC1,G0C1,GC2,G0C2,GC3,G0C3,GC4,G0C4,GC5,G0C5;	// contacts
				GCstamp(G,0,1,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptC_T,GC0,G0C0);
				GCstamp(G,3,4,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptC_T,GC1,G0C1);	
				GCstamp(G,6,7,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptC_T,GC2,G0C2);
				GCstamp(G,9,10,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptC_T,GC3,G0C3);	
				GCstamp(G,12,13,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptC_T,GC4,G0C4);
				GCstamp(G,15,16,output_magnet_factor * OutputMagnetLength - output_magnet_deduction,ioptM_W,ioptC_T,GC5,G0C5);

				Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3,GF4,G0F4,GF5,G0F5;	// magnets
				GFstamp(G,1,2,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptM_T,GF0,G0F0);
				GFstamp(G,4,5,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptM_T,GF1,G0F1);	
				GFstamp(G,7,8,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptM_T,GF2,G0F2);
				GFstamp(G,10,11,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptM_T,GF3,G0F3);	
				GFstamp(G,13,14,input_magnet_factor * InputMagnetLength - input_magnet_deduction,ioptM_W,ioptM_T,GF4,G0F4);
				GFstamp(G,16,17,output_magnet_factor * OutputMagnetLength - output_magnet_deduction,ioptM_W,ioptM_T,GF5,G0F5);

				Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3,GN4,G0N4,GN5,G0N5;	// channels
				GNstamp(G,2,18,midC_W,midC_T,midC_L,GN0,G0N0);
				GNstamp(G,5,18,midC_W,midC_T,midC_L,GN1,G0N1);
				GNstamp(G,8,18,midC_W,midC_T,midC_L,GN2,G0N2);
				GNstamp(G,11,18,midC_W,midC_T,midC_L,GN3,G0N3);
				GNstamp(G,14,18,midC_W,midC_T,midC_L,GN4,G0N4);
				GNstamp(G,17,18,midC_W,midC_T,midC_L,GN5,G0N5);
				// **** for now, we use global parameters to define ground dimensions for simplicity of codes ****
				Matrix2d GN6,G0N6,GN7,G0N7,GN8,G0N8,GN9,G0N9,GN10,G0N10;	// grounds
				GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN6,G0N6);
				GGstamp(G,5,-1,iptG_W,iptG_T,iptG_L,GN7,G0N7);
				GGstamp(G,8,-1,iptG_W,iptG_T,iptG_L,GN8,G0N8);
				GGstamp(G,11,-1,iptG_W,iptG_T,iptG_L,GN9,G0N9);
				GGstamp(G,14,-1,iptG_W,iptG_T,iptG_L,GN10,G0N10);

				Vstamp(G,0,19);
				Vstamp(G,3,20);
				Vstamp(G,6,21);
				Vstamp(G,9,22);
				Vstamp(G,12,23);
				Vstamp(G,15,24);

				VectorXd b(2*(n+6));
				b << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				     Vdd, 0, -Vdd, 0, Vdd, 0, -Vdd, 0, Vdd, 0, 0, 0;

				VectorXd x_V = G.inverse() * b;	
			//	cout<<G<<endl;
			//	cout<<x_V<<endl;
				MatrixXd r_V = ReshapeV(x_V);	
			//	cout<<r_V<<endl;	

				Vector2d I1, I7, I13, I16;
				I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
				I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

				I7(0) = GF2(0,0)*(r_V(0,7)-r_V(0,8)) + GF2(0,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(0,0)*r_V(0,8) - G0F2(0,1)*r_V(1,8);
				I7(1) = GF2(1,0)*(r_V(0,7)-r_V(0,8)) + GF2(1,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(1,0)*r_V(0,8) - G0F2(1,1)*r_V(1,8);

				I13(0) = GF4(0,0)*(r_V(0,13)-r_V(0,14)) + GF4(0,1)*(r_V(1,13)-r_V(1,14)); // - G0F4(0,0)*r_V(0,13) - G0F4(0,1)*r_V(1,14);
				I13(1) = GF4(1,0)*(r_V(0,13)-r_V(0,14)) + GF4(1,1)*(r_V(1,13)-r_V(1,14)); // - G0F4(1,0)*r_V(0,13) - G0F4(1,1)*r_V(1,14);
			
				I16(0) = GF5(0,0)*(r_V(0,17)-r_V(0,16)) + GF5(0,1)*(r_V(1,17)-r_V(1,16)); // - G0F5(0,0)*r_V(0,16) - G0F5(0,1)*r_V(1,16);
				I16(1) = GF5(1,0)*(r_V(0,17)-r_V(0,16)) + GF5(1,1)*(r_V(1,17)-r_V(1,16)); // - G0F5(1,0)*r_V(0,16) - G0F5(1,1)*r_V(1,16);			

			//	cout<<I10(1)<<endl;
				
				double SIE = I16(1)/I1(0);
			//	cout <<SIE<<endl;

				double R_C = rho_C*ioptC_T/(0.5*ioptM_W*InputMagnetLength);	// contact resistance considered
				double R_F = rho_F*ioptM_T/(0.5*ioptM_W*InputMagnetLength);	// ferromagnet resistance considered
				double R_N = rho_N*midC_T/(0.5*ioptM_W*InputMagnetLength);	// channel resistance considered
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);	// ground resistance considered
				double Ic_real = Vdd/( R_C + R_F + R_N + R_G );

				double Ns = Ms*OutputMagnetLength*ioptM_W*ioptM_T/mu_B;
				t_sw5to1 = 2*fsw*q*Ns/abs(Ic_real*SIE);

			}
			else if (swoption == 1) 
			{		
				t_sw5to1 = 1000e-9;

				cout<<"Error: switching time in the InternalDelayPreCalc function for 5to1 of swoption==1 not properly defined!!!"<<endl;
				assert(false);
			}
			else if (swoption == 0) 
			{
				int n=14;	// # of nodes
				MatrixXd G = MatrixXd::Zero(2*n,2*n);	// initialize with zeros every time, otherwise it seems it will accumulate??

				Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3,GF4,G0F4,GF5,G0F5;
				GFstamp(G,0,1,InputMagnetLength,ioptM_W,ioptM_T,GF0,G0F0);
				GFstamp(G,2,3,InputMagnetLength,ioptM_W,ioptM_T,GF1,G0F1);	// set the size of fixed magnet same as the first input magnet
				GFstamp(G,4,5,InputMagnetLength,ioptM_W,ioptM_T,GF2,G0F2);
				GFstamp(G,6,7,InputMagnetLength,ioptM_W,ioptM_T,GF3,G0F3);
				GFstamp(G,8,9,InputMagnetLength,ioptM_W,ioptM_T,GF4,G0F4);
				GFstamp(G,10,11,OutputMagnetLength,ioptM_W,ioptM_T,GF5,G0F5);
			
				Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3,GN4,G0N4,GN5,G0N5;
				GNstamp(G,1,13,midC_W,midC_T,midC_L,GN0,G0N0);
				GNstamp(G,3,13,midC_W,midC_T,midC_L,GN1,G0N1); // set the size of input chanenl of fixed magnet same as the channel of first input magnet
				GNstamp(G,5,13,midC_W,midC_T,midC_L,GN2,G0N2);
				GNstamp(G,7,13,midC_W,midC_T,midC_L,GN3,G0N3);
				GNstamp(G,9,13,midC_W,midC_T,midC_L,GN4,G0N4);
				GNstamp(G,11,13,midC_W,midC_T,midC_L,GN5,G0N5); // the channel length after joint set as 50e-9
		
				Matrix2d GN6,G0N6,GN7,G0N7,GN8,G0N8,GN9,G0N9,GN10,G0N10,GN11,G0N11;
				GNstamp(G,1,-1,midC_W,midC_T,100*lambda_N,GN6,G0N6);
				GNstamp(G,3,-1,midC_W,midC_T,100*lambda_N,GN7,G0N7);
				GNstamp(G,5,-1,midC_W,midC_T,100*lambda_N,GN8,G0N8);
				GNstamp(G,7,-1,midC_W,midC_T,100*lambda_N,GN9,G0N9); 
				GNstamp(G,9,-1,midC_W,midC_T,100*lambda_N,GN10,G0N10);
				GNstamp(G,11,12,midC_W,midC_T,100*lambda_N,GN11,G0N11);

				VectorXd C(2*n);
				C << Ic, 0, 0, 0, -Ic, 0, 0, 0, Ic, 0, 0, 0, -Ic, 0, 0, 0, Ic, 0, 0, 0,// worse case delay should be considered, therefore two positve Ic and one negative Ic is used here
				       0, 0, 0, 0, 0, 0, 0, 0;

				VectorXd x_V = G.inverse() * C;	
			//	cout<<G<<endl;
			//	cout<<x_V<<endl;
				MatrixXd r_V = ReshapeV(x_V);	
			//	cout<<r_V<<endl;	
			
				Vector2d I0, I4, I8, I10;
				I0(0) = GF0(0,0)*(r_V(0,0)-r_V(0,1)) + GF0(0,1)*(r_V(1,0)-r_V(1,1)) - G0F0(0,0)*r_V(0,1) - G0F0(0,1)*r_V(1,1);
				I0(1) = GF0(1,0)*(r_V(0,0)-r_V(0,1)) + GF0(1,1)*(r_V(1,0)-r_V(1,1)) - G0F0(1,0)*r_V(0,1) - G0F0(1,1)*r_V(1,1);

				I4(0) = GF2(0,0)*(r_V(0,4)-r_V(0,5)) + GF2(0,1)*(r_V(1,4)-r_V(1,5)) - G0F2(0,0)*r_V(0,5) - G0F2(0,1)*r_V(1,5);
				I4(1) = GF2(1,0)*(r_V(0,4)-r_V(0,5)) + GF2(1,1)*(r_V(1,4)-r_V(1,5)) - G0F2(1,0)*r_V(0,5) - G0F2(1,1)*r_V(1,5);
			
				I8(0) = GF4(0,0)*(r_V(0,8)-r_V(0,9)) + GF4(0,1)*(r_V(1,8)-r_V(1,9)) - G0F4(0,0)*r_V(0,9) - G0F4(0,1)*r_V(1,9);
				I8(1) = GF4(1,0)*(r_V(0,8)-r_V(0,9)) + GF4(1,1)*(r_V(1,8)-r_V(1,9)) - G0F4(1,0)*r_V(0,9) - G0F4(1,1)*r_V(1,9);

				I10(0) = GF5(0,0)*(r_V(0,10)-r_V(0,11)) + GF5(0,1)*(r_V(1,10)-r_V(1,11)) - G0F5(0,0)*r_V(0,11) - G0F5(0,1)*r_V(1,11);
				I10(1) = GF5(1,0)*(r_V(0,10)-r_V(0,11)) + GF5(1,1)*(r_V(1,10)-r_V(1,11)) - G0F5(1,0)*r_V(0,11) - G0F5(1,1)*r_V(1,11);

				double SIE = I10(1)/I0(0);
				// half of the area for the magnet and channel under the source is considered
				double Ic_real = Vdd/(rho_F*ioptM_T/(0.5*ioptM_W*InputMagnetLength) + rho_N*midC_T/(0.5*ioptM_W*InputMagnetLength));
		
				double Ns = Ms*OutputMagnetLength*ioptM_W*ioptM_T/mu_B;
				t_sw5to1 = 2*fsw*q*Ns/abs(Ic_real*SIE);

			}

			Delay3to1Table[CurrentPair] = t_sw3to1;
			Delay5to1Table[CurrentPair] = t_sw5to1;

		//	cout<<"Current Input Length & Output Length: ("<<CurrentPair.first<<","<<CurrentPair.second<<"):"<<Delay3to1Table[CurrentPair]<<" and "<<Delay5to1Table[CurrentPair]<<endl;
			j += 1;
		} 
		i += 1;
	}

}

void IntDelayCalcMatch(gate *igate_ptr, map<pair<int, int>, double> &Delay3to1Table, map<pair<int, int>, double> &Delay5to1Table)
{
	// calculate the delay inside the gate
	if ( std::find(gatetype3to1list.begin(),gatetype3to1list.end(),igate_ptr->gatetype) != gatetype3to1list.end() )
	{
		pair<int, int> CurrentPair ((int)(igate_ptr->iptM_L[0]*1e9), (int)(igate_ptr->optM_L*1e9));
		double t_sw = Delay3to1Table[CurrentPair];
		
		igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
		fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);
	//	if (igate_ptr->name == "n9") cout<<"n9: ("<<igate_ptr->iptM_L[0]<<", "<<igate_ptr->optM_L<<"): delay "<<t_sw<<endl;

	}
	else if ( std::find(gatetype5to1list.begin(),gatetype5to1list.end(),igate_ptr->gatetype) != gatetype5to1list.end() )
	{
		pair<int, int> CurrentPair ((int)(igate_ptr->iptM_L[0]*1e9), (int)(igate_ptr->optM_L*1e9));
		double t_sw = Delay5to1Table[CurrentPair];
		
		igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
		fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

	}
	else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),igate_ptr->gatetype) != gatetype1to1list.end() )
	{

/* The following part is commented because for Invertor/Buffer (not matter they are in the netlist initially or they are added through BufferInsertion process.
Their delay should be set as zero inside the gate. And external delay should be calculated only by InterconnectDelayCalc() using their own ogate and previous
gates' ogate dimensions. */
		double t_sw = 0;

		igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
		fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

	}
	else if ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),igate_ptr->gatetype) != gatetypepriinputlist.end() )
	{
		double t_sw = 0;

		igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
		fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

	}
	else if ( std::find(gatetypeconstlist.begin(),gatetypeconstlist.end(),igate_ptr->gatetype) != gatetypeconstlist.end() )
	{
		double t_sw = 0;
		
		igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
		fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

	}
	else 
	{	
		cout<<"Not Defined GateType: "<<igate_ptr->gatetype<<endl;
		assert(false);
	}
}

void IntDelayCalc(gate *igate_ptr)
{
	// calculate the delay inside the gate
	if ( std::find(gatetype3to1list.begin(),gatetype3to1list.end(),igate_ptr->gatetype) != gatetype3to1list.end() )
	{
		if (swoption == 3)	// modeling 2 includes contacts on every magnet and grounds under input side 
		{	// rules to write the netlist same as in swoption 0 or 1

			int n=13;
			MatrixXd G = MatrixXd::Zero(2*(n+4),2*(n+4));
			// **** for now, we use the magnet dimensions and global parameters to define contact dimensions for simplicity of codes ****
			Matrix2d GC0,G0C0,GC1,G0C1,GC2,G0C2,GC3,G0C3;	// contacts
			GCstamp(G,0,1,input_contact_factor * igate_ptr->iptM_L[0] + input_contact_adjustment,igate_ptr->iptM_W[0],ioptC_T,GC0,G0C0);
			GCstamp(G,3,4,input_contact_factor * igate_ptr->iptM_L[0] + input_contact_adjustment,igate_ptr->iptM_W[0],ioptC_T,GC1,G0C1);	
			GCstamp(G,6,7,input_contact_factor * igate_ptr->iptM_L[1] + input_contact_adjustment,igate_ptr->iptM_W[1],ioptC_T,GC2,G0C2);
			GCstamp(G,9,10,output_contact_factor * igate_ptr->optM_L + output_contact_adjustment,igate_ptr->optM_W,ioptC_T,GC3,G0C3);

			Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3;	// magnets
			GFstamp(G,1,2,input_magnet_factor * igate_ptr->iptM_L[0] + input_magnet_adjustment,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF0,G0F0);
			GFstamp(G,4,5,input_magnet_factor * igate_ptr->iptM_L[0] + input_magnet_adjustment,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF1,G0F1);	
			GFstamp(G,7,8,input_magnet_factor * igate_ptr->iptM_L[1] + input_magnet_adjustment,igate_ptr->iptM_W[1],igate_ptr->iptM_T[1],GF2,G0F2);
			GFstamp(G,10,11,output_magnet_factor * igate_ptr->optM_L + output_magnet_adjustment,igate_ptr->optM_W,igate_ptr->optM_T,GF3,G0F3);

			Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3;	// channels
			GNstamp(G,2,12,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN0,G0N0);
			GNstamp(G,5,12,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN1,G0N1);
			GNstamp(G,8,12,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],igate_ptr->iptC_L[1],GN2,G0N2);
			GNstamp(G,11,12,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[1],GN3,G0N3);
			// **** for now, we use global parameters to define ground dimensions for simplicity of codes ****
			Matrix2d GN4,G0N4,GN5,G0N5,GN6,G0N6;	// grounds
			GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN4,G0N4);
			GGstamp(G,5,-1,iptG_W,iptG_T,iptG_L,GN5,G0N5);
			GGstamp(G,8,-1,iptG_W,iptG_T,iptG_L,GN6,G0N6);			

			Vstamp(G,0,13);
			Vstamp(G,3,14);
			Vstamp(G,6,15);
			Vstamp(G,9,16);

			VectorXd b(2*(n+4));
			b << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     Vdd, 0, -Vdd, 0, Vdd, 0, 0, 0;

			VectorXd x_V = G.inverse() * b;	
		//	cout<<G<<endl;
		//	cout<<x_V<<endl;
			MatrixXd r_V = ReshapeV(x_V);	
		//	cout<<r_V<<endl;	

			Vector2d I1, I7, I10;
			I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
			I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

			I7(0) = GF2(0,0)*(r_V(0,7)-r_V(0,8)) + GF2(0,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(0,0)*r_V(0,8) - G0F2(0,1)*r_V(1,8);
			I7(1) = GF2(1,0)*(r_V(0,7)-r_V(0,8)) + GF2(1,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(1,0)*r_V(0,8) - G0F2(1,1)*r_V(1,8);
			
			I10(0) = GF3(0,0)*(r_V(0,11)-r_V(0,10)) + GF3(0,1)*(r_V(1,11)-r_V(1,10)); // - G0F3(0,0)*r_V(0,10) - G0F3(0,1)*r_V(1,10);
			I10(1) = GF3(1,0)*(r_V(0,11)-r_V(0,10)) + GF3(1,1)*(r_V(1,11)-r_V(1,10)); // - G0F3(1,0)*r_V(0,10) - G0F3(1,1)*r_V(1,10);			

		//	cout<<I10(1)<<endl;
				
			double SIE = I10(1)/I1(0);
		//	cout <<SIE<<endl;

			double R_C = rho_C*ioptC_T/(igate_ptr->iptM_W[0]*(input_contact_factor * igate_ptr->iptM_L[0] + input_contact_adjustment));	// contact resistance considered
			double R_F = rho_F*igate_ptr->iptM_T[0]/(igate_ptr->iptM_W[0]*(input_magnet_factor * igate_ptr->iptM_L[0] + input_magnet_adjustment));	// ferromagnet resistance considered
			double R_N = rho_N*igate_ptr->iptC_T[0]/(igate_ptr->iptM_W[0]*(input_nonmagnet_factor * igate_ptr->iptM_L[0] + input_nonmagnet_adjustment));	// channel resistance considered
			double R_G = rho_G*iptG_L/(iptG_W*iptG_T);	// ground resistance considered
			double Ic_real = Vdd/( R_C + R_F + R_N + R_G );

			double Ns = Ms*igate_ptr->optM_L*igate_ptr->optM_W*igate_ptr->optM_T/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);

			igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
			fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

		}
		else if (swoption == 2)	// modeling 2 includes contacts on every magnet and grounds under input side 
		{	// rules to write the netlist same as in swoption 0 or 1

			int n=13;
			MatrixXd G = MatrixXd::Zero(2*(n+4),2*(n+4));
			// **** for now, we use the magnet dimensions and global parameters to define contact dimensions for simplicity of codes ****
			Matrix2d GC0,G0C0,GC1,G0C1,GC2,G0C2,GC3,G0C3;	// contacts
			GCstamp(G,0,1,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],ioptC_T,GC0,G0C0);
			GCstamp(G,3,4,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],ioptC_T,GC1,G0C1);	
			GCstamp(G,6,7,input_magnet_factor * igate_ptr->iptM_L[1] - input_magnet_deduction,igate_ptr->iptM_W[1],ioptC_T,GC2,G0C2);
			GCstamp(G,9,10,output_magnet_factor * igate_ptr->optM_L - output_magnet_deduction,igate_ptr->optM_W,ioptC_T,GC3,G0C3);

			Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3;	// magnets
			GFstamp(G,1,2,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF0,G0F0);
			GFstamp(G,4,5,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF1,G0F1);	
			GFstamp(G,7,8,input_magnet_factor * igate_ptr->iptM_L[1] - input_magnet_deduction,igate_ptr->iptM_W[1],igate_ptr->iptM_T[1],GF2,G0F2);
			GFstamp(G,10,11,output_magnet_factor * igate_ptr->optM_L - output_magnet_deduction,igate_ptr->optM_W,igate_ptr->optM_T,GF3,G0F3);

			Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3;	// channels
			GNstamp(G,2,12,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN0,G0N0);
			GNstamp(G,5,12,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN1,G0N1);
			GNstamp(G,8,12,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],igate_ptr->iptC_L[1],GN2,G0N2);
			GNstamp(G,11,12,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[1],GN3,G0N3);
			// **** for now, we use global parameters to define ground dimensions for simplicity of codes ****
			Matrix2d GN4,G0N4,GN5,G0N5,GN6,G0N6;	// grounds
			GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN4,G0N4);
			GGstamp(G,5,-1,iptG_W,iptG_T,iptG_L,GN5,G0N5);
			GGstamp(G,8,-1,iptG_W,iptG_T,iptG_L,GN6,G0N6);			

			Vstamp(G,0,13);
			Vstamp(G,3,14);
			Vstamp(G,6,15);
			Vstamp(G,9,16);

			VectorXd b(2*(n+4));
			b << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     Vdd, 0, -Vdd, 0, Vdd, 0, 0, 0;

			VectorXd x_V = G.inverse() * b;	
		//	cout<<G<<endl;
		//	cout<<x_V<<endl;
			MatrixXd r_V = ReshapeV(x_V);	
		//	cout<<r_V<<endl;	

			Vector2d I1, I7, I10;
			I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
			I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

			I7(0) = GF2(0,0)*(r_V(0,7)-r_V(0,8)) + GF2(0,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(0,0)*r_V(0,8) - G0F2(0,1)*r_V(1,8);
			I7(1) = GF2(1,0)*(r_V(0,7)-r_V(0,8)) + GF2(1,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(1,0)*r_V(0,8) - G0F2(1,1)*r_V(1,8);
			
			I10(0) = GF3(0,0)*(r_V(0,11)-r_V(0,10)) + GF3(0,1)*(r_V(1,11)-r_V(1,10)); // - G0F3(0,0)*r_V(0,10) - G0F3(0,1)*r_V(1,10);
			I10(1) = GF3(1,0)*(r_V(0,11)-r_V(0,10)) + GF3(1,1)*(r_V(1,11)-r_V(1,10)); // - G0F3(1,0)*r_V(0,10) - G0F3(1,1)*r_V(1,10);			

		//	cout<<I10(1)<<endl;
				
			double SIE = I10(1)/I1(0);
		//	cout <<SIE<<endl;

			double R_C = rho_C*ioptC_T/(0.5*igate_ptr->iptM_W[0]*igate_ptr->iptM_L[0]);	// contact resistance considered
			double R_F = rho_F*igate_ptr->iptM_T[0]/(0.5*igate_ptr->iptM_W[0]*igate_ptr->iptM_L[0]);	// ferromagnet resistance considered
			double R_N = rho_N*igate_ptr->iptC_T[0]/(0.5*igate_ptr->iptM_W[0]*igate_ptr->iptM_L[0]);	// channel resistance considered
			double R_G = rho_G*iptG_L/(iptG_W*iptG_T);	// ground resistance considered
			double Ic_real = Vdd/( R_C + R_F + R_N + R_G );

			double Ns = Ms*igate_ptr->optM_L*igate_ptr->optM_W*igate_ptr->optM_T/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);

			igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
			fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

		}
		else if (swoption == 1) 
		{	// depending on we choose to impose constant current source or voltage source, the matrices will be different

			int n=10;	// # of nodes
			MatrixXd G = MatrixXd::Zero(2*(n+3),2*(n+3));	// initialize with zeros every time, otherwise it seems it will accumulate??

			Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3;
			GFstamp(G,0,1,igate_ptr->iptM_L[0],igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF0,G0F0);
			GFstamp(G,2,3,igate_ptr->iptM_L[0],igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF1,G0F1);	// set the size of fixed magnet same as the first input magnet
			GFstamp(G,4,5,igate_ptr->iptM_L[1],igate_ptr->iptM_W[1],igate_ptr->iptM_T[1],GF2,G0F2);
			GFstamp(G,6,7,igate_ptr->optM_L,igate_ptr->optM_W,igate_ptr->optM_T,GF3,G0F3);
			
			Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3;
			GNstamp(G,1,8,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN0,G0N0);
			GNstamp(G,3,8,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN1,G0N1);	// set the size of input chanenl of fixed magnet same as the channel of first input magnet
			GNstamp(G,5,8,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],igate_ptr->iptC_L[1],GN2,G0N2);
			GNstamp(G,7,8,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[1],GN3,G0N3);	// the channel length after joint set as 50e-9
			
			Matrix2d GN4,G0N4,GN5,G0N5,GN6,G0N6,GN7,G0N7;
			GNstamp(G,1,-1,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],100*lambda_N,GN4,G0N4);
			GNstamp(G,3,-1,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],100*lambda_N,GN5,G0N5);
			GNstamp(G,5,-1,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],100*lambda_N,GN6,G0N6);
			GNstamp(G,7,9,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],100*lambda_N,GN7,G0N7); 

			Vstamp(G,0,10);	// two input numbers, first is the ID/order of the node attached to the voltage source, second is how many nodes already there in the G matrix
			Vstamp(G,2,11);
			Vstamp(G,4,12);

			VectorXd b(2*(n+3));
			b << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     Vdd, 0, -Vdd, 0, Vdd, 0;

			VectorXd x_V = G.inverse() * b;	
		//	cout<<G<<endl;
		//	cout<<x_V<<endl;
			MatrixXd r_V = ReshapeV(x_V);	
		//	cout<<r_V<<endl;	

			
			Vector2d I0, I4, I6;
			I0(0) = GF0(0,0)*(r_V(0,0)-r_V(0,1)) + GF0(0,1)*(r_V(1,0)-r_V(1,1)) - G0F0(0,0)*r_V(0,1) - G0F0(0,1)*r_V(1,1);
			I0(1) = GF0(1,0)*(r_V(0,0)-r_V(0,1)) + GF0(1,1)*(r_V(1,0)-r_V(1,1)) - G0F0(1,0)*r_V(0,1) - G0F0(1,1)*r_V(1,1);

			I4(0) = GF2(0,0)*(r_V(0,4)-r_V(0,5)) + GF2(0,1)*(r_V(1,4)-r_V(1,5)) - G0F2(0,0)*r_V(0,5) - G0F2(0,1)*r_V(1,5);
			I4(1) = GF2(1,0)*(r_V(0,4)-r_V(0,5)) + GF2(1,1)*(r_V(1,4)-r_V(1,5)) - G0F2(1,0)*r_V(0,5) - G0F2(1,1)*r_V(1,5);
			
			I6(0) = GF3(0,0)*(r_V(0,6)-r_V(0,7)) + GF3(0,1)*(r_V(1,6)-r_V(1,7)) - G0F3(0,0)*r_V(0,7) - G0F3(0,1)*r_V(1,7);
			I6(1) = GF3(1,0)*(r_V(0,6)-r_V(0,7)) + GF3(1,1)*(r_V(1,6)-r_V(1,7)) - G0F3(1,0)*r_V(0,7) - G0F3(1,1)*r_V(1,7);
		
		//	cout<<I6(1)<<endl;
				
			double SIE = I6(1)/I0(0);
		//	cout <<SIE<<endl;
			
			double Ns = Ms*igate_ptr->optM_L*igate_ptr->optM_W*igate_ptr->optM_T/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(I6(1));
			
			igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
			fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

		}
		else if (swoption == 0) 
		{
		
			int n=10;	// # of nodes
			MatrixXd G = MatrixXd::Zero(2*n,2*n);	// initialize with zeros every time, otherwise it seems it will accumulate??

			Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3;
			GFstamp(G,0,1,igate_ptr->iptM_L[0],igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF0,G0F0);
			GFstamp(G,2,3,igate_ptr->iptM_L[0],igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF1,G0F1);	// set the size of fixed magnet same as the first input magnet
			GFstamp(G,4,5,igate_ptr->iptM_L[1],igate_ptr->iptM_W[1],igate_ptr->iptM_T[1],GF2,G0F2);
			GFstamp(G,6,7,igate_ptr->optM_L,igate_ptr->optM_W,igate_ptr->optM_T,GF3,G0F3);
			
			Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3;
			GNstamp(G,1,8,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN0,G0N0);
			GNstamp(G,3,8,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN1,G0N1);	// set the size of input chanenl of fixed magnet same as the channel of first input magnet
			GNstamp(G,5,8,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],igate_ptr->iptC_L[1],GN2,G0N2);
			GNstamp(G,7,8,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[1],GN3,G0N3);	// the channel length after joint set as 50e-9
			
			Matrix2d GN4,G0N4,GN5,G0N5,GN6,G0N6,GN7,G0N7;
			GNstamp(G,1,-1,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],100*lambda_N,GN4,G0N4);
			GNstamp(G,3,-1,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],100*lambda_N,GN5,G0N5);
			GNstamp(G,5,-1,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],100*lambda_N,GN6,G0N6);
			GNstamp(G,7,9,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],100*lambda_N,GN7,G0N7); 

			VectorXd C(2*n);
			C << Ic, 0, 0, 0, -Ic, 0, 0, 0, Ic, 0, 0, 0, 	// worse case delay should be considered, therefore two positve Ic and one negative Ic is used here
			       0, 0, 0, 0, 0, 0, 0, 0;

			VectorXd x_V = G.inverse() * C;	
		//	cout<<G<<endl;
		//	cout<<x_V<<endl;
			MatrixXd r_V = ReshapeV(x_V);	
		//	cout<<r_V<<endl;	

			Vector2d I0, I4, I6;
			I0(0) = GF0(0,0)*(r_V(0,0)-r_V(0,1)) + GF0(0,1)*(r_V(1,0)-r_V(1,1)) - G0F0(0,0)*r_V(0,1) - G0F0(0,1)*r_V(1,1);
			I0(1) = GF0(1,0)*(r_V(0,0)-r_V(0,1)) + GF0(1,1)*(r_V(1,0)-r_V(1,1)) - G0F0(1,0)*r_V(0,1) - G0F0(1,1)*r_V(1,1);

			I4(0) = GF2(0,0)*(r_V(0,4)-r_V(0,5)) + GF2(0,1)*(r_V(1,4)-r_V(1,5)) - G0F2(0,0)*r_V(0,5) - G0F2(0,1)*r_V(1,5);
			I4(1) = GF2(1,0)*(r_V(0,4)-r_V(0,5)) + GF2(1,1)*(r_V(1,4)-r_V(1,5)) - G0F2(1,0)*r_V(0,5) - G0F2(1,1)*r_V(1,5);
			
			I6(0) = GF3(0,0)*(r_V(0,6)-r_V(0,7)) + GF3(0,1)*(r_V(1,6)-r_V(1,7)) - G0F3(0,0)*r_V(0,7) - G0F3(0,1)*r_V(1,7);
			I6(1) = GF3(1,0)*(r_V(0,6)-r_V(0,7)) + GF3(1,1)*(r_V(1,6)-r_V(1,7)) - G0F3(1,0)*r_V(0,7) - G0F3(1,1)*r_V(1,7);
		
			double SIE = I6(1)/I0(0);
			// half of the area for the magnet and channel under the source is considered		
			double Ic_real = Vdd/(rho_F*igate_ptr->iptM_T[0]/(0.5*igate_ptr->iptM_W[0]*igate_ptr->iptM_L[0]) + rho_N*igate_ptr->iptC_T[0]/(0.5*igate_ptr->iptM_W[0]*igate_ptr->iptM_L[0]));

			double Ns = Ms*igate_ptr->optM_L*igate_ptr->optM_W*igate_ptr->optM_T/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);
			
			igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
			fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

		}
	}
	else if ( std::find(gatetype5to1list.begin(),gatetype5to1list.end(),igate_ptr->gatetype) != gatetype5to1list.end() )
	{
		if (swoption == 3)
		{
			int n=19;
			MatrixXd G = MatrixXd::Zero(2*(n+6),2*(n+6));
			// **** for now, we use the magnet dimensions and global parameters to define contact dimensions for simplicity of codes ****
			Matrix2d GC0,G0C0,GC1,G0C1,GC2,G0C2,GC3,G0C3,GC4,G0C4,GC5,G0C5;	// contacts
			GCstamp(G,0,1,input_contact_factor * igate_ptr->iptM_L[0] + input_contact_adjustment,igate_ptr->iptM_W[0],ioptC_T,GC0,G0C0);
			GCstamp(G,3,4,input_contact_factor * igate_ptr->iptM_L[0] + input_contact_adjustment,igate_ptr->iptM_W[0],ioptC_T,GC1,G0C1);	
			GCstamp(G,6,7,input_contact_factor * igate_ptr->iptM_L[0] + input_contact_adjustment,igate_ptr->iptM_W[0],ioptC_T,GC2,G0C2);
			GCstamp(G,9,10,input_contact_factor * igate_ptr->iptM_L[0] + input_contact_adjustment,igate_ptr->iptM_W[0],ioptC_T,GC3,G0C3);	
			GCstamp(G,12,13,input_contact_factor * igate_ptr->iptM_L[1] + input_contact_adjustment,igate_ptr->iptM_W[1],ioptC_T,GC4,G0C4);
			GCstamp(G,15,16,output_contact_factor * igate_ptr->optM_L + output_contact_adjustment,igate_ptr->optM_W,ioptC_T,GC5,G0C5);

			Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3,GF4,G0F4,GF5,G0F5;	// magnets
			GFstamp(G,1,2,input_magnet_factor * igate_ptr->iptM_L[0] + input_magnet_adjustment,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF0,G0F0);
			GFstamp(G,4,5,input_magnet_factor * igate_ptr->iptM_L[0] + input_magnet_adjustment,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF1,G0F1);	
			GFstamp(G,7,8,input_magnet_factor * igate_ptr->iptM_L[0] + input_magnet_adjustment,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF2,G0F2);
			GFstamp(G,10,11,input_magnet_factor * igate_ptr->iptM_L[0] + input_magnet_adjustment,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF3,G0F3);	
			GFstamp(G,13,14,input_magnet_factor * igate_ptr->iptM_L[1] + input_magnet_adjustment,igate_ptr->iptM_W[1],igate_ptr->iptM_T[1],GF4,G0F4);
			GFstamp(G,16,17,output_magnet_factor * igate_ptr->optM_L + output_magnet_adjustment,igate_ptr->optM_W,igate_ptr->optM_T,GF5,G0F5);

			Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3,GN4,G0N4,GN5,G0N5;	// channels
			GNstamp(G,2,18,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN0,G0N0);
			GNstamp(G,5,18,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN1,G0N1);
			GNstamp(G,8,18,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN2,G0N2);
			GNstamp(G,11,18,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN3,G0N3);
			GNstamp(G,14,18,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],igate_ptr->iptC_L[1],GN4,G0N4);
			GNstamp(G,17,18,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[1],GN5,G0N5);
			// **** for now, we use global parameters to define ground dimensions for simplicity of codes ****
			Matrix2d GN6,G0N6,GN7,G0N7,GN8,G0N8,GN9,G0N9,GN10,G0N10;	// grounds
			GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN6,G0N6);
			GGstamp(G,5,-1,iptG_W,iptG_T,iptG_L,GN7,G0N7);
			GGstamp(G,8,-1,iptG_W,iptG_T,iptG_L,GN8,G0N8);
			GGstamp(G,11,-1,iptG_W,iptG_T,iptG_L,GN9,G0N9);
			GGstamp(G,14,-1,iptG_W,iptG_T,iptG_L,GN10,G0N10);

			Vstamp(G,0,19);
			Vstamp(G,3,20);
			Vstamp(G,6,21);
			Vstamp(G,9,22);
			Vstamp(G,12,23);
			Vstamp(G,15,24);

			VectorXd b(2*(n+6));
			b << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     Vdd, 0, -Vdd, 0, Vdd, 0, -Vdd, 0, Vdd, 0, 0, 0;

			VectorXd x_V = G.inverse() * b;	
		//	cout<<G<<endl;
		//	cout<<x_V<<endl;
			MatrixXd r_V = ReshapeV(x_V);	
		//	cout<<r_V<<endl;	

			Vector2d I1, I7, I13, I16;
			I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
			I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

			I7(0) = GF2(0,0)*(r_V(0,7)-r_V(0,8)) + GF2(0,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(0,0)*r_V(0,8) - G0F2(0,1)*r_V(1,8);
			I7(1) = GF2(1,0)*(r_V(0,7)-r_V(0,8)) + GF2(1,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(1,0)*r_V(0,8) - G0F2(1,1)*r_V(1,8);

			I13(0) = GF4(0,0)*(r_V(0,13)-r_V(0,14)) + GF4(0,1)*(r_V(1,13)-r_V(1,14)); // - G0F4(0,0)*r_V(0,13) - G0F4(0,1)*r_V(1,14);
			I13(1) = GF4(1,0)*(r_V(0,13)-r_V(0,14)) + GF4(1,1)*(r_V(1,13)-r_V(1,14)); // - G0F4(1,0)*r_V(0,13) - G0F4(1,1)*r_V(1,14);
			
			I16(0) = GF5(0,0)*(r_V(0,17)-r_V(0,16)) + GF5(0,1)*(r_V(1,17)-r_V(1,16)); // - G0F5(0,0)*r_V(0,16) - G0F5(0,1)*r_V(1,16);
			I16(1) = GF5(1,0)*(r_V(0,17)-r_V(0,16)) + GF5(1,1)*(r_V(1,17)-r_V(1,16)); // - G0F5(1,0)*r_V(0,16) - G0F5(1,1)*r_V(1,16);			

		//	cout<<I10(1)<<endl;
				
			double SIE = I16(1)/I1(0);
		//	cout <<SIE<<endl;

			double R_C = rho_C*ioptC_T/(igate_ptr->iptM_W[0]*(input_contact_factor * igate_ptr->iptM_L[0] + input_contact_adjustment));	// contact resistance considered
			double R_F = rho_F*igate_ptr->iptM_T[0]/(igate_ptr->iptM_W[0]*(input_magnet_factor * igate_ptr->iptM_L[0] + input_magnet_adjustment));	// ferromagnet resistance considered
			double R_N = rho_N*igate_ptr->iptC_T[0]/(igate_ptr->iptM_W[0]*(input_nonmagnet_factor * igate_ptr->iptM_L[0] + input_nonmagnet_adjustment));	// channel resistance considered
			double R_G = rho_G*iptG_L/(iptG_W*iptG_T);	// ground resistance considered
			double Ic_real = Vdd/( R_C + R_F + R_N + R_G );

			double Ns = Ms*igate_ptr->optM_L*igate_ptr->optM_W*igate_ptr->optM_T/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);

			igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
			fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);
		
		}
		else if (swoption == 2)
		{
			int n=19;
			MatrixXd G = MatrixXd::Zero(2*(n+6),2*(n+6));
			// **** for now, we use the magnet dimensions and global parameters to define contact dimensions for simplicity of codes ****
			Matrix2d GC0,G0C0,GC1,G0C1,GC2,G0C2,GC3,G0C3,GC4,G0C4,GC5,G0C5;	// contacts
			GCstamp(G,0,1,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],ioptC_T,GC0,G0C0);
			GCstamp(G,3,4,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],ioptC_T,GC1,G0C1);	
			GCstamp(G,6,7,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],ioptC_T,GC2,G0C2);
			GCstamp(G,9,10,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],ioptC_T,GC3,G0C3);	
			GCstamp(G,12,13,input_magnet_factor * igate_ptr->iptM_L[1] - input_magnet_deduction,igate_ptr->iptM_W[1],ioptC_T,GC4,G0C4);
			GCstamp(G,15,16,output_magnet_factor * igate_ptr->optM_L - output_magnet_deduction,igate_ptr->optM_W,ioptC_T,GC5,G0C5);

			Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3,GF4,G0F4,GF5,G0F5;	// magnets
			GFstamp(G,1,2,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF0,G0F0);
			GFstamp(G,4,5,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF1,G0F1);	
			GFstamp(G,7,8,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF2,G0F2);
			GFstamp(G,10,11,input_magnet_factor * igate_ptr->iptM_L[0] - input_magnet_deduction,igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF3,G0F3);	
			GFstamp(G,13,14,input_magnet_factor * igate_ptr->iptM_L[1] - input_magnet_deduction,igate_ptr->iptM_W[1],igate_ptr->iptM_T[1],GF4,G0F4);
			GFstamp(G,16,17,output_magnet_factor * igate_ptr->optM_L - output_magnet_deduction,igate_ptr->optM_W,igate_ptr->optM_T,GF5,G0F5);

			Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3,GN4,G0N4,GN5,G0N5;	// channels
			GNstamp(G,2,18,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN0,G0N0);
			GNstamp(G,5,18,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN1,G0N1);
			GNstamp(G,8,18,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN2,G0N2);
			GNstamp(G,11,18,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN3,G0N3);
			GNstamp(G,14,18,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],igate_ptr->iptC_L[1],GN4,G0N4);
			GNstamp(G,17,18,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[1],GN5,G0N5);
			// **** for now, we use global parameters to define ground dimensions for simplicity of codes ****
			Matrix2d GN6,G0N6,GN7,G0N7,GN8,G0N8,GN9,G0N9,GN10,G0N10;	// grounds
			GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN6,G0N6);
			GGstamp(G,5,-1,iptG_W,iptG_T,iptG_L,GN7,G0N7);
			GGstamp(G,8,-1,iptG_W,iptG_T,iptG_L,GN8,G0N8);
			GGstamp(G,11,-1,iptG_W,iptG_T,iptG_L,GN9,G0N9);
			GGstamp(G,14,-1,iptG_W,iptG_T,iptG_L,GN10,G0N10);

			Vstamp(G,0,19);
			Vstamp(G,3,20);
			Vstamp(G,6,21);
			Vstamp(G,9,22);
			Vstamp(G,12,23);
			Vstamp(G,15,24);

			VectorXd b(2*(n+6));
			b << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			     Vdd, 0, -Vdd, 0, Vdd, 0, -Vdd, 0, Vdd, 0, 0, 0;

			VectorXd x_V = G.inverse() * b;	
		//	cout<<G<<endl;
		//	cout<<x_V<<endl;
			MatrixXd r_V = ReshapeV(x_V);	
		//	cout<<r_V<<endl;	

			Vector2d I1, I7, I13, I16;
			I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
			I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

			I7(0) = GF2(0,0)*(r_V(0,7)-r_V(0,8)) + GF2(0,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(0,0)*r_V(0,8) - G0F2(0,1)*r_V(1,8);
			I7(1) = GF2(1,0)*(r_V(0,7)-r_V(0,8)) + GF2(1,1)*(r_V(1,7)-r_V(1,8)); // - G0F2(1,0)*r_V(0,8) - G0F2(1,1)*r_V(1,8);

			I13(0) = GF4(0,0)*(r_V(0,13)-r_V(0,14)) + GF4(0,1)*(r_V(1,13)-r_V(1,14)); // - G0F4(0,0)*r_V(0,13) - G0F4(0,1)*r_V(1,14);
			I13(1) = GF4(1,0)*(r_V(0,13)-r_V(0,14)) + GF4(1,1)*(r_V(1,13)-r_V(1,14)); // - G0F4(1,0)*r_V(0,13) - G0F4(1,1)*r_V(1,14);
			
			I16(0) = GF5(0,0)*(r_V(0,17)-r_V(0,16)) + GF5(0,1)*(r_V(1,17)-r_V(1,16)); // - G0F5(0,0)*r_V(0,16) - G0F5(0,1)*r_V(1,16);
			I16(1) = GF5(1,0)*(r_V(0,17)-r_V(0,16)) + GF5(1,1)*(r_V(1,17)-r_V(1,16)); // - G0F5(1,0)*r_V(0,16) - G0F5(1,1)*r_V(1,16);			

		//	cout<<I10(1)<<endl;
				
			double SIE = I16(1)/I1(0);
		//	cout <<SIE<<endl;

			double R_C = rho_C*ioptC_T/(0.5*igate_ptr->iptM_W[0]*igate_ptr->iptM_L[0]);	// contact resistance considered
			double R_F = rho_F*igate_ptr->iptM_T[0]/(0.5*igate_ptr->iptM_W[0]*igate_ptr->iptM_L[0]);	// ferromagnet resistance considered
			double R_N = rho_N*igate_ptr->iptC_T[0]/(0.5*igate_ptr->iptM_W[0]*igate_ptr->iptM_L[0]);	// channel resistance considered
			double R_G = rho_G*iptG_L/(iptG_W*iptG_T);	// ground resistance considered
			double Ic_real = Vdd/( R_C + R_F + R_N + R_G );

			double Ns = Ms*igate_ptr->optM_L*igate_ptr->optM_W*igate_ptr->optM_T/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);

			igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
			fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);
		
		}
		else if (swoption == 1) 
		{		
			double t_sw = 1000e-9;

			igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
			fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

		}
		else if (swoption == 0) 
		{
			int n=14;	// # of nodes
			MatrixXd G = MatrixXd::Zero(2*n,2*n);	// initialize with zeros every time, otherwise it seems it will accumulate??

			Matrix2d GF0,G0F0,GF1,G0F1,GF2,G0F2,GF3,G0F3,GF4,G0F4,GF5,G0F5;
			GFstamp(G,0,1,igate_ptr->iptM_L[0],igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF0,G0F0);
			GFstamp(G,2,3,igate_ptr->iptM_L[0],igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF1,G0F1);	// set the size of fixed magnet same as the first input magnet
			GFstamp(G,4,5,igate_ptr->iptM_L[1],igate_ptr->iptM_W[1],igate_ptr->iptM_T[1],GF2,G0F2);
			GFstamp(G,6,7,igate_ptr->iptM_L[0],igate_ptr->iptM_W[0],igate_ptr->iptM_T[0],GF3,G0F3);
			GFstamp(G,8,9,igate_ptr->iptM_L[1],igate_ptr->iptM_W[1],igate_ptr->iptM_T[1],GF4,G0F4);
			GFstamp(G,10,11,igate_ptr->optM_L,igate_ptr->optM_W,igate_ptr->optM_T,GF5,G0F5);
			
			Matrix2d GN0,G0N0,GN1,G0N1,GN2,G0N2,GN3,G0N3,GN4,G0N4,GN5,G0N5;
			GNstamp(G,1,13,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN0,G0N0);
			GNstamp(G,3,13,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN1,G0N1); // set the size of input chanenl of fixed magnet same as the channel of first input magnet
			GNstamp(G,5,13,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],igate_ptr->iptC_L[1],GN2,G0N2);
			GNstamp(G,7,13,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],igate_ptr->iptC_L[0],GN3,G0N3);
			GNstamp(G,9,13,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],igate_ptr->iptC_L[1],GN4,G0N4);
			GNstamp(G,11,13,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],igate_ptr->iptC_L[1],GN5,G0N5); // the channel length after joint set as 50e-9
		
			Matrix2d GN6,G0N6,GN7,G0N7,GN8,G0N8,GN9,G0N9,GN10,G0N10,GN11,G0N11;
			GNstamp(G,1,-1,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],100*lambda_N,GN6,G0N6);
			GNstamp(G,3,-1,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],100*lambda_N,GN7,G0N7);
			GNstamp(G,5,-1,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],100*lambda_N,GN8,G0N8);
			GNstamp(G,7,-1,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],100*lambda_N,GN9,G0N9); 
			GNstamp(G,9,-1,igate_ptr->iptC_W[1],igate_ptr->iptC_T[1],100*lambda_N,GN10,G0N10);
			GNstamp(G,11,12,igate_ptr->iptC_W[0],igate_ptr->iptC_T[0],100*lambda_N,GN11,G0N11);

			VectorXd C(2*n);
			C << Ic, 0, 0, 0, -Ic, 0, 0, 0, Ic, 0, 0, 0, -Ic, 0, 0, 0, Ic, 0, 0, 0,// worse case delay should be considered, therefore two positve Ic and one negative Ic is used here
			       0, 0, 0, 0, 0, 0, 0, 0;

			VectorXd x_V = G.inverse() * C;	
		//	cout<<G<<endl;
		//	cout<<x_V<<endl;
			MatrixXd r_V = ReshapeV(x_V);	
		//	cout<<r_V<<endl;	
			
			Vector2d I0, I4, I8, I10;
			I0(0) = GF0(0,0)*(r_V(0,0)-r_V(0,1)) + GF0(0,1)*(r_V(1,0)-r_V(1,1)) - G0F0(0,0)*r_V(0,1) - G0F0(0,1)*r_V(1,1);
			I0(1) = GF0(1,0)*(r_V(0,0)-r_V(0,1)) + GF0(1,1)*(r_V(1,0)-r_V(1,1)) - G0F0(1,0)*r_V(0,1) - G0F0(1,1)*r_V(1,1);

			I4(0) = GF2(0,0)*(r_V(0,4)-r_V(0,5)) + GF2(0,1)*(r_V(1,4)-r_V(1,5)) - G0F2(0,0)*r_V(0,5) - G0F2(0,1)*r_V(1,5);
			I4(1) = GF2(1,0)*(r_V(0,4)-r_V(0,5)) + GF2(1,1)*(r_V(1,4)-r_V(1,5)) - G0F2(1,0)*r_V(0,5) - G0F2(1,1)*r_V(1,5);
			
			I8(0) = GF4(0,0)*(r_V(0,8)-r_V(0,9)) + GF4(0,1)*(r_V(1,8)-r_V(1,9)) - G0F4(0,0)*r_V(0,9) - G0F4(0,1)*r_V(1,9);
			I8(1) = GF4(1,0)*(r_V(0,8)-r_V(0,9)) + GF4(1,1)*(r_V(1,8)-r_V(1,9)) - G0F4(1,0)*r_V(0,9) - G0F4(1,1)*r_V(1,9);

			I10(0) = GF5(0,0)*(r_V(0,10)-r_V(0,11)) + GF5(0,1)*(r_V(1,10)-r_V(1,11)) - G0F5(0,0)*r_V(0,11) - G0F5(0,1)*r_V(1,11);
			I10(1) = GF5(1,0)*(r_V(0,10)-r_V(0,11)) + GF5(1,1)*(r_V(1,10)-r_V(1,11)) - G0F5(1,0)*r_V(0,11) - G0F5(1,1)*r_V(1,11);

			double SIE = I10(1)/I0(0);
			// half of the area for the magnet and channel under the source is considered
			double Ic_real = Vdd/(rho_F*igate_ptr->iptM_T[0]/(0.5*igate_ptr->iptM_W[0]*igate_ptr->iptM_L[0]) + rho_N*igate_ptr->iptC_T[0]/(0.5*igate_ptr->iptM_W[0]*igate_ptr->iptM_L[0]));
		
			double Ns = Ms*igate_ptr->optM_L*igate_ptr->optM_W*igate_ptr->optM_T/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);
			
			igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
			fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);
		}
	}
	else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),igate_ptr->gatetype) != gatetype1to1list.end() )
	{

/* The following part is commented because for Invertor/Buffer (not matter they are in the netlist initially or they are added through BufferInsertion process.
Their delay should be set as zero inside the gate. And external delay should be calculated only by InterconnectDelayCalc() using their own ogate and previous
gates' ogate dimensions. */
		double t_sw = 0;

		igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
		fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

	}
	else if ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),igate_ptr->gatetype) != gatetypepriinputlist.end() )
	{
		double t_sw = 0;

		igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
		fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

	}
	else if ( std::find(gatetypeconstlist.begin(),gatetypeconstlist.end(),igate_ptr->gatetype) != gatetypeconstlist.end() )
	{
		double t_sw = 0;
		
		igate_ptr->int_ipt_delay.resize(igate_ptr->ipt_nodes.size());
		fill(igate_ptr->int_ipt_delay.begin(), igate_ptr->int_ipt_delay.end(), t_sw);

	}
	else 
	{	
		cout<<"Not Defined GateType: "<<igate_ptr->gatetype<<endl;
		assert(false);
	}
	
}

void InterconnectDelayCalc(gate *igate_ptr, gate *ogate_ptr, int ipt_order)
{
/* identify the output gatetype ogate first, if it is buffer (not matter inserted or original ones), then the output magnet size for 
the output gate is used, otherwise for the 3/5 input gates, their corresponding input magnet size is used. Other type or gates will 
only give error */
	if ( ( std::find(gatetype3to1list.begin(),gatetype3to1list.end(),ogate_ptr->gatetype) != gatetype3to1list.end() ) || ( std::find(gatetype5to1list.begin(),gatetype5to1list.end(),ogate_ptr->gatetype) != gatetype5to1list.end() ) ) 
	{
		if (swoption == 3)
		{
			int n = 6;
			MatrixXd G = MatrixXd::Zero(2*(n+2),2*(n+2));
			// **** for now, we use the magnet dimensions and global parameters to define contact dimensions for simplicity of codes ****
			Matrix2d GC0,G0C0,GC1,G0C1;	// contacts
			GCstamp(G,0,1,input_contact_factor * igate_ptr->optM_L + input_contact_adjustment,igate_ptr->optM_W,ioptC_T,GC0,G0C0);
			GCstamp(G,3,4,output_contact_factor * ogate_ptr->iptM_L[ipt_order] + output_contact_adjustment,ogate_ptr->iptM_W[ipt_order],ioptC_T,GC1,G0C1);			

			Matrix2d GF0,G0F0,GF1,G0F1;
			GFstamp(G,1,2,input_magnet_factor * igate_ptr->optM_L + input_magnet_adjustment,igate_ptr->optM_W,igate_ptr->optM_T,GF0,G0F0);
			GFstamp(G,4,5,output_magnet_factor * ogate_ptr->iptM_L[ipt_order] + output_magnet_adjustment,ogate_ptr->iptM_W[ipt_order],ogate_ptr->iptM_T[ipt_order],GF1,G0F1);

			Matrix2d GN0,G0N0;
			GNstamp(G,2,5,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],ogate_ptr->ipt_dis[ipt_order]*1e-9,GN0,G0N0);

			Matrix2d GN1,G0N1;
			GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN1,G0N1);

			Vstamp(G,0,6);
			Vstamp(G,3,7);
		//	cout<<G<<endl;

			VectorXd b(2*(n+2));

			b<<0, 0, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 0,
			   Vdd, 0, 0, 0;

			VectorXd x_V = G.inverse() * b;

			MatrixXd r_V = ReshapeV(x_V);
		//	cout<<r_V<<endl;

			Vector2d I1, I4;
			I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
			I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

			I4(0) = GF1(0,0)*(r_V(0,5)-r_V(0,4)) + GF1(0,1)*(r_V(1,5)-r_V(1,4)); // - G0F1(0,0)*r_V(0,4) - G0F1(0,1)*r_V(1,4);
			I4(1) = GF1(1,0)*(r_V(0,5)-r_V(0,4)) + GF1(1,1)*(r_V(1,5)-r_V(1,4)); // - G0F1(1,0)*r_V(0,4) - G0F1(1,1)*r_V(1,4);

			double SIE = abs(I4(1)/I1(0));
			double Ic_real; // depending on the input is primary input or not, the current calculation should be different in considering 0.5*Area or not
			if (igate_ptr->gatetype != "PRI_INPUT")
			{
				double R_C = rho_C*ioptC_T/(igate_ptr->optM_W*(input_contact_factor * igate_ptr->optM_L + input_contact_adjustment));
				double R_F = rho_F*igate_ptr->optM_T/(igate_ptr->optM_W*(input_magnet_factor * igate_ptr->optM_L + input_magnet_adjustment));
				double R_N = rho_N*ogate_ptr->iptC_T[ipt_order]/(igate_ptr->optM_W*(input_nonmagnet_factor * igate_ptr->optM_L + input_nonmagnet_adjustment));
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_real = Vdd/( R_C + R_F + R_N + R_G );
			}
			else if (igate_ptr->gatetype == "PRI_INPUT")
			{
				double R_C = rho_C*ioptC_T/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_F = rho_F*igate_ptr->optM_T/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_N = rho_N*ogate_ptr->iptC_T[ipt_order]/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_real = Vdd/( R_C + R_F + R_N + R_G );
			}
		
			double Ns = Ms*ogate_ptr->iptM_L[ipt_order]*ogate_ptr->iptM_W[ipt_order]*ogate_ptr->iptM_T[ipt_order]/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);
		//	cout<<t_sw<<endl;
		
			ogate_ptr->ext_ipt_delay[ipt_order] = t_sw;

		}
		else if (swoption == 2)
		{
			int n = 6;
			MatrixXd G = MatrixXd::Zero(2*(n+2),2*(n+2));
			// **** for now, we use the magnet dimensions and global parameters to define contact dimensions for simplicity of codes ****
			Matrix2d GC0,G0C0,GC1,G0C1;	// contacts
			GCstamp(G,0,1,input_magnet_factor * igate_ptr->optM_L - input_magnet_deduction,igate_ptr->optM_W,ioptC_T,GC0,G0C0);
			GCstamp(G,3,4,output_magnet_factor * ogate_ptr->iptM_L[ipt_order] - output_magnet_deduction,ogate_ptr->iptM_W[ipt_order],ioptC_T,GC1,G0C1);			

			Matrix2d GF0,G0F0,GF1,G0F1;
			GFstamp(G,1,2,input_magnet_factor * igate_ptr->optM_L - input_magnet_deduction,igate_ptr->optM_W,igate_ptr->optM_T,GF0,G0F0);
			GFstamp(G,4,5,output_magnet_factor * ogate_ptr->iptM_L[ipt_order] - output_magnet_deduction,ogate_ptr->iptM_W[ipt_order],ogate_ptr->iptM_T[ipt_order],GF1,G0F1);

			Matrix2d GN0,G0N0;
			GNstamp(G,2,5,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],ogate_ptr->ipt_dis[ipt_order]*1e-9,GN0,G0N0);

			Matrix2d GN1,G0N1;
			GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN1,G0N1);

			Vstamp(G,0,6);
			Vstamp(G,3,7);
		//	cout<<G<<endl;

			VectorXd b(2*(n+2));

			b<<0, 0, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 0,
			   Vdd, 0, 0, 0;

			VectorXd x_V = G.inverse() * b;

			MatrixXd r_V = ReshapeV(x_V);
		//	cout<<r_V<<endl;

			Vector2d I1, I4;
			I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
			I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

			I4(0) = GF1(0,0)*(r_V(0,5)-r_V(0,4)) + GF1(0,1)*(r_V(1,5)-r_V(1,4)); // - G0F1(0,0)*r_V(0,4) - G0F1(0,1)*r_V(1,4);
			I4(1) = GF1(1,0)*(r_V(0,5)-r_V(0,4)) + GF1(1,1)*(r_V(1,5)-r_V(1,4)); // - G0F1(1,0)*r_V(0,4) - G0F1(1,1)*r_V(1,4);

			double SIE = abs(I4(1)/I1(0));
			double Ic_real; // depending on the input is primary input or not, the current calculation should be different in considering 0.5*Area or not
			if (igate_ptr->gatetype != "PRI_INPUT")
			{
				double R_C = rho_C*ioptC_T/(0.5*igate_ptr->optM_W*igate_ptr->optM_L);
				double R_F = rho_F*igate_ptr->optM_T/(0.5*igate_ptr->optM_W*igate_ptr->optM_L);
				double R_N = rho_N*ogate_ptr->iptC_T[ipt_order]/(0.5*igate_ptr->optM_W*igate_ptr->optM_L);
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_real = Vdd/( R_C + R_F + R_N + R_G );
			}
			else if (igate_ptr->gatetype == "PRI_INPUT")
			{
				double R_C = rho_C*ioptC_T/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_F = rho_F*igate_ptr->optM_T/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_N = rho_N*ogate_ptr->iptC_T[ipt_order]/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_real = Vdd/( R_C + R_F + R_N + R_G );
			}
		
			double Ns = Ms*ogate_ptr->iptM_L[ipt_order]*ogate_ptr->iptM_W[ipt_order]*ogate_ptr->iptM_T[ipt_order]/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);
		//	cout<<t_sw<<endl;
		
			ogate_ptr->ext_ipt_delay[ipt_order] = t_sw;

		}
		else if (swoption == 1) 
		{
			int n = 5;
			MatrixXd G = MatrixXd::Zero(2*(n+1),2*(n+1));

			Matrix2d GF0,G0F0,GF1,G0F1;
			GFstamp(G,0,1,igate_ptr->optM_L,igate_ptr->optM_W,igate_ptr->optM_T,GF0,G0F0);
			GFstamp(G,2,3,ogate_ptr->iptM_L[ipt_order],ogate_ptr->iptM_W[ipt_order],ogate_ptr->iptM_T[ipt_order],GF1,G0F1);

			Matrix2d GN0,G0N0;
			GNstamp(G,1,3,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],ogate_ptr->ipt_dis[ipt_order]*1e-9,GN0,G0N0);
		//	cout<<ogate_ptr->iptC_W[ipt_order]<<endl;
		//	cout<<ogate_ptr->iptC_T[ipt_order]<<endl;
		//	cout<<ogate_ptr->ipt_dis[ipt_order]<<endl;

			Matrix2d GN1,G0N1,GN2,G0N2;
			GNstamp(G,1,-1,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],100*lambda_N,GN1,G0N1);
			GNstamp(G,3,4,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],100*lambda_N,GN2,G0N2);

			Vstamp(G,0,5);
		//	cout<<G<<endl;

			VectorXd b(2*(n+1));

			b<<0, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 
			   Vdd, 0;

			VectorXd x_V = G.inverse() * b;

			MatrixXd r_V = ReshapeV(x_V);
		//	cout<<r_V<<endl;

			Vector2d I0, I3;
			I0(0) = GF0(0,0)*(r_V(0,0)-r_V(0,1)) + GF0(0,1)*(r_V(1,0)-r_V(1,1)) - G0F0(0,0)*r_V(0,1) - G0F0(0,1)*r_V(1,1);
			I0(1) = GF0(1,0)*(r_V(0,0)-r_V(0,1)) + GF0(1,1)*(r_V(1,0)-r_V(1,1)) - G0F0(1,0)*r_V(0,1) - G0F0(1,1)*r_V(1,1);

			I3(0) = GF1(0,0)*(r_V(0,2)-r_V(0,3)) + GF1(0,1)*(r_V(1,2)-r_V(1,3)) - G0F1(0,0)*r_V(0,3) - G0F1(0,1)*r_V(1,3);
			I3(1) = GF1(1,0)*(r_V(0,2)-r_V(0,3)) + GF1(1,1)*(r_V(1,2)-r_V(1,3)) - G0F1(1,0)*r_V(0,3) - G0F1(1,1)*r_V(1,3);

			double SIE = abs(I3(1)/I0(0));
		//	cout<<"SIE: "<<SIE<<endl;
				
			double Ns = Ms*ogate_ptr->iptM_L[ipt_order]*ogate_ptr->iptM_W[ipt_order]*ogate_ptr->iptM_T[ipt_order]/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(I3(1));
		//	cout<<t_sw<<endl;
		
			ogate_ptr->ext_ipt_delay[ipt_order] = t_sw;
		//	cout<<ogate_ptr->ext_ipt_delay[ipt_order]<<endl;

		}
		else if (swoption == 0) 
		{
			int n = 5;
			MatrixXd G = MatrixXd::Zero(2*n,2*n);

			Matrix2d GF0,G0F0,GF1,G0F1;
			GFstamp(G,0,1,igate_ptr->optM_L,igate_ptr->optM_W,igate_ptr->optM_T,GF0,G0F0);
			GFstamp(G,2,3,ogate_ptr->iptM_L[ipt_order],ogate_ptr->iptM_W[ipt_order],ogate_ptr->iptM_T[ipt_order],GF1,G0F1);

			Matrix2d GN0,G0N0;
			GNstamp(G,1,3,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],ogate_ptr->ipt_dis[ipt_order]*1e-9,GN0,G0N0);

			Matrix2d GN1,G0N1,GN2,G0N2;
			GNstamp(G,1,-1,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],100*lambda_N,GN1,G0N1);
			GNstamp(G,3,4,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],100*lambda_N,GN2,G0N2);

			VectorXd C(2*n);

			C<<Ic, 0, 0, 0, 0, 0, 0, 0, 0, 0;

			VectorXd x_V = G.inverse() * C;

			MatrixXd r_V = ReshapeV(x_V);
		//	cout<<r_V<<endl;

			Vector2d I0, I3;
			I0(0) = GF0(0,0)*(r_V(0,0)-r_V(0,1)) + GF0(0,1)*(r_V(1,0)-r_V(1,1)) - G0F0(0,0)*r_V(0,1) - G0F0(0,1)*r_V(1,1);
			I0(1) = GF0(1,0)*(r_V(0,0)-r_V(0,1)) + GF0(1,1)*(r_V(1,0)-r_V(1,1)) - G0F0(1,0)*r_V(0,1) - G0F0(1,1)*r_V(1,1);

			I3(0) = GF1(0,0)*(r_V(0,2)-r_V(0,3)) + GF1(0,1)*(r_V(1,2)-r_V(1,3)) - G0F1(0,0)*r_V(0,3) - G0F1(0,1)*r_V(1,3);
			I3(1) = GF1(1,0)*(r_V(0,2)-r_V(0,3)) + GF1(1,1)*(r_V(1,2)-r_V(1,3)) - G0F1(1,0)*r_V(0,3) - G0F1(1,1)*r_V(1,3);

			double SIE = abs(I3(1)/I0(0));
			double Ic_real; // depending on the input is primary input or not, the current calculation should be different in considering 0.5*Area or not
			if (igate_ptr->gatetype != "PRI_INPUT")
			{	
				Ic_real = Vdd/(rho_F*igate_ptr->optM_T/(0.5*igate_ptr->optM_W*igate_ptr->optM_L) + rho_N*ogate_ptr->iptC_T[ipt_order]/(0.5*igate_ptr->optM_W*igate_ptr->optM_L));
			}
			else if (igate_ptr->gatetype == "PRI_INPUT")
			{	
				Ic_real = Vdd/(rho_F*igate_ptr->optM_T/(igate_ptr->optM_W*igate_ptr->optM_L) + rho_N*ogate_ptr->iptC_T[ipt_order]/(igate_ptr->optM_W*igate_ptr->optM_L));
			}
		
			double Ns = Ms*ogate_ptr->iptM_L[ipt_order]*ogate_ptr->iptM_W[ipt_order]*ogate_ptr->iptM_T[ipt_order]/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);
		//	cout<<t_sw<<endl;
		
			ogate_ptr->ext_ipt_delay[ipt_order] = t_sw;

		}

	}
	else if ( std::find(gatetype1to1list.begin(),gatetype1to1list.end(),ogate_ptr->gatetype) != gatetype1to1list.end() )
	{
		if (swoption == 3)
		{
			int n = 6;
			MatrixXd G = MatrixXd::Zero(2*(n+2),2*(n+2));

			Matrix2d GC0,G0C0,GC1,G0C1;
			GCstamp(G,0,1,input_contact_factor * igate_ptr->optM_L + input_contact_adjustment,igate_ptr->optM_W,ioptC_T,GC0,G0C0);
			GCstamp(G,3,4,output_contact_factor * ogate_ptr->optM_L + output_contact_adjustment,ogate_ptr->optM_W,ioptC_T,GC1,G0C1);

			Matrix2d GF0,G0F0,GF1,G0F1;
			GFstamp(G,1,2,input_magnet_factor * igate_ptr->optM_L + input_magnet_adjustment,igate_ptr->optM_W,igate_ptr->optM_T,GF0,G0F0);
			GFstamp(G,4,5,output_magnet_factor * ogate_ptr->optM_L + output_magnet_adjustment,ogate_ptr->optM_W,ogate_ptr->optM_T,GF1,G0F1);

			Matrix2d GN0,G0N0;
			GNstamp(G,2,5,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],ogate_ptr->ipt_dis[ipt_order]*1e-9,GN0,G0N0);

			Matrix2d GN1,G0N1;
			GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN1,G0N1);

			Vstamp(G,0,6);
			Vstamp(G,3,7);
		//	cout<<G<<endl;

			VectorXd b(2*(n+2));

			b<<0, 0, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 0,
			   Vdd, 0, 0, 0;

			VectorXd x_V = G.inverse() * b;

			MatrixXd r_V = ReshapeV(x_V);
		//	cout<<r_V<<endl;

			Vector2d I1, I4;
			I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
			I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

			I4(0) = GF1(0,0)*(r_V(0,5)-r_V(0,4)) + GF1(0,1)*(r_V(1,5)-r_V(1,4)); // - G0F1(0,0)*r_V(0,4) - G0F1(0,1)*r_V(1,4);
			I4(1) = GF1(1,0)*(r_V(0,5)-r_V(0,4)) + GF1(1,1)*(r_V(1,5)-r_V(1,4)); // - G0F1(1,0)*r_V(0,4) - G0F1(1,1)*r_V(1,4);

			double SIE = abs(I4(1)/I1(0));
		//	cout<<"SIE: "<<SIE<<endl;
			double Ic_real;
			if (igate_ptr->gatetype != "PRI_INPUT")
			{	
				double R_C = rho_C*ioptC_T/(igate_ptr->optM_W*(input_contact_factor * igate_ptr->optM_L + input_contact_adjustment));
				double R_F = rho_F*igate_ptr->optM_T/(igate_ptr->optM_W*(input_magnet_factor * igate_ptr->optM_L + input_magnet_adjustment));
				double R_N = rho_N*ogate_ptr->iptC_T[ipt_order]/(igate_ptr->optM_W*(input_nonmagnet_factor * igate_ptr->optM_L + input_nonmagnet_adjustment));
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_real = Vdd/( R_C + R_F + R_N + R_G );
			}
			else if (igate_ptr->gatetype == "PRI_INPUT")
			{
				double R_C = rho_C*ioptC_T/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_F = rho_F*igate_ptr->optM_T/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_N = rho_N*ogate_ptr->iptC_T[ipt_order]/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_real = Vdd/( R_C + R_F + R_N + R_G );
			}
	
			double Ns = Ms*ogate_ptr->optM_L*ogate_ptr->optM_W*ogate_ptr->optM_T/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);
		//	cout<<t_sw<<endl;
		
			ogate_ptr->ext_ipt_delay[ipt_order] = t_sw;
		//	cout<<ogate_ptr->ext_ipt_delay[ipt_order]<<endl;

		}
		else if (swoption == 2)
		{
			int n = 6;
			MatrixXd G = MatrixXd::Zero(2*(n+2),2*(n+2));

			Matrix2d GC0,G0C0,GC1,G0C1;
			GCstamp(G,0,1,input_magnet_factor * igate_ptr->optM_L - input_magnet_deduction,igate_ptr->optM_W,ioptC_T,GC0,G0C0);
			GCstamp(G,3,4,output_magnet_factor * ogate_ptr->optM_L - output_magnet_deduction,ogate_ptr->optM_W,ioptC_T,GC1,G0C1);

			Matrix2d GF0,G0F0,GF1,G0F1;
			GFstamp(G,1,2,input_magnet_factor * igate_ptr->optM_L - input_magnet_deduction,igate_ptr->optM_W,igate_ptr->optM_T,GF0,G0F0);
			GFstamp(G,4,5,output_magnet_factor * ogate_ptr->optM_L - output_magnet_deduction,ogate_ptr->optM_W,ogate_ptr->optM_T,GF1,G0F1);

			Matrix2d GN0,G0N0;
			GNstamp(G,2,5,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],ogate_ptr->ipt_dis[ipt_order]*1e-9,GN0,G0N0);

			Matrix2d GN1,G0N1;
			GGstamp(G,2,-1,iptG_W,iptG_T,iptG_L,GN1,G0N1);

			Vstamp(G,0,6);
			Vstamp(G,3,7);
		//	cout<<G<<endl;

			VectorXd b(2*(n+2));

			b<<0, 0, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 0,
			   Vdd, 0, 0, 0;

			VectorXd x_V = G.inverse() * b;

			MatrixXd r_V = ReshapeV(x_V);
		//	cout<<r_V<<endl;

			Vector2d I1, I4;
			I1(0) = GF0(0,0)*(r_V(0,1)-r_V(0,2)) + GF0(0,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(0,0)*r_V(0,2) - G0F0(0,1)*r_V(1,2);
			I1(1) = GF0(1,0)*(r_V(0,1)-r_V(0,2)) + GF0(1,1)*(r_V(1,1)-r_V(1,2)); // - G0F0(1,0)*r_V(0,2) - G0F0(1,1)*r_V(1,2);

			I4(0) = GF1(0,0)*(r_V(0,5)-r_V(0,4)) + GF1(0,1)*(r_V(1,5)-r_V(1,4)); // - G0F1(0,0)*r_V(0,4) - G0F1(0,1)*r_V(1,4);
			I4(1) = GF1(1,0)*(r_V(0,5)-r_V(0,4)) + GF1(1,1)*(r_V(1,5)-r_V(1,4)); // - G0F1(1,0)*r_V(0,4) - G0F1(1,1)*r_V(1,4);

			double SIE = abs(I4(1)/I1(0));
		//	cout<<"SIE: "<<SIE<<endl;
			double Ic_real;
			if (igate_ptr->gatetype != "PRI_INPUT")
			{	
				double R_C = rho_C*ioptC_T/(0.5*igate_ptr->optM_W*igate_ptr->optM_L);
				double R_F = rho_F*igate_ptr->optM_T/(0.5*igate_ptr->optM_W*igate_ptr->optM_L);
				double R_N = rho_N*ogate_ptr->iptC_T[ipt_order]/(0.5*igate_ptr->optM_W*igate_ptr->optM_L);
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_real = Vdd/( R_C + R_F + R_N + R_G );
			}
			else if (igate_ptr->gatetype == "PRI_INPUT")
			{
				double R_C = rho_C*ioptC_T/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_F = rho_F*igate_ptr->optM_T/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_N = rho_N*ogate_ptr->iptC_T[ipt_order]/(igate_ptr->optM_W*igate_ptr->optM_L);
				double R_G = rho_G*iptG_L/(iptG_W*iptG_T);
				Ic_real = Vdd/( R_C + R_F + R_N + R_G );
			}
	
			double Ns = Ms*ogate_ptr->optM_L*ogate_ptr->optM_W*ogate_ptr->optM_T/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);
		//	cout<<t_sw<<endl;
		
			ogate_ptr->ext_ipt_delay[ipt_order] = t_sw;
		//	cout<<ogate_ptr->ext_ipt_delay[ipt_order]<<endl;

		}
		else if (swoption == 1) 
		{
			int n = 5;
			MatrixXd G = MatrixXd::Zero(2*(n+1),2*(n+1));

			Matrix2d GF0,G0F0,GF1,G0F1;
			GFstamp(G,0,1,igate_ptr->optM_L,igate_ptr->optM_W,igate_ptr->optM_T,GF0,G0F0);
			GFstamp(G,2,3,ogate_ptr->optM_L,ogate_ptr->optM_W,ogate_ptr->optM_T,GF1,G0F1);

			Matrix2d GN0,G0N0;
			GNstamp(G,1,3,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],ogate_ptr->ipt_dis[ipt_order]*1e-9,GN0,G0N0);
		//	cout<<ogate_ptr->iptC_W[ipt_order]<<endl;
		//	cout<<ogate_ptr->iptC_T[ipt_order]<<endl;
		//	cout<<ogate_ptr->ipt_dis[ipt_order]<<endl;

			Matrix2d GN1,G0N1,GN2,G0N2;
			GNstamp(G,1,-1,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],100*lambda_N,GN1,G0N1);
			GNstamp(G,3,4,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],100*lambda_N,GN2,G0N2);

			Vstamp(G,0,5);
		//	cout<<G<<endl;

			VectorXd b(2*(n+1));

			b<<0, 0, 0, 0, 0,
			   0, 0, 0, 0, 0, 
			   Vdd, 0;

			VectorXd x_V = G.inverse() * b;

			MatrixXd r_V = ReshapeV(x_V);
		//	cout<<r_V<<endl;

			Vector2d I0, I3;
			I0(0) = GF0(0,0)*(r_V(0,0)-r_V(0,1)) + GF0(0,1)*(r_V(1,0)-r_V(1,1)) - G0F0(0,0)*r_V(0,1) - G0F0(0,1)*r_V(1,1);
			I0(1) = GF0(1,0)*(r_V(0,0)-r_V(0,1)) + GF0(1,1)*(r_V(1,0)-r_V(1,1)) - G0F0(1,0)*r_V(0,1) - G0F0(1,1)*r_V(1,1);

			I3(0) = GF1(0,0)*(r_V(0,2)-r_V(0,3)) + GF1(0,1)*(r_V(1,2)-r_V(1,3)) - G0F1(0,0)*r_V(0,3) - G0F1(0,1)*r_V(1,3);
			I3(1) = GF1(1,0)*(r_V(0,2)-r_V(0,3)) + GF1(1,1)*(r_V(1,2)-r_V(1,3)) - G0F1(1,0)*r_V(0,3) - G0F1(1,1)*r_V(1,3);

			double SIE = abs(I3(1)/I0(0));
		//	cout<<"SIE: "<<SIE<<endl;
				
			double Ns = Ms*ogate_ptr->optM_L*ogate_ptr->optM_W*ogate_ptr->optM_T/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(I3(1));
		//	cout<<t_sw<<endl;
		
			ogate_ptr->ext_ipt_delay[ipt_order] = t_sw;
		//	cout<<ogate_ptr->ext_ipt_delay[ipt_order]<<endl;

		}
		else if (swoption == 0) 
		{
			int n = 5;
			MatrixXd G = MatrixXd::Zero(2*n,2*n);

			Matrix2d GF0,G0F0,GF1,G0F1;
			GFstamp(G,0,1,igate_ptr->optM_L,igate_ptr->optM_W,igate_ptr->optM_T,GF0,G0F0);
			GFstamp(G,2,3,ogate_ptr->optM_L,ogate_ptr->optM_W,ogate_ptr->optM_T,GF1,G0F1);

			Matrix2d GN0,G0N0;
			GNstamp(G,1,3,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],ogate_ptr->ipt_dis[ipt_order]*1e-9,GN0,G0N0);

			Matrix2d GN1,G0N1,GN2,G0N2;
			GNstamp(G,1,-1,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],100*lambda_N,GN1,G0N1);
			GNstamp(G,3,4,ogate_ptr->iptC_W[ipt_order],ogate_ptr->iptC_T[ipt_order],100*lambda_N,GN2,G0N2);

			VectorXd C(2*n);

			C<<Ic, 0, 0, 0, 0, 0, 0, 0, 0, 0;

			VectorXd x_V = G.inverse() * C;

			MatrixXd r_V = ReshapeV(x_V);
		//	cout<<r_V<<endl;

			Vector2d I0, I3;
			I0(0) = GF0(0,0)*(r_V(0,0)-r_V(0,1)) + GF0(0,1)*(r_V(1,0)-r_V(1,1)) - G0F0(0,0)*r_V(0,1) - G0F0(0,1)*r_V(1,1);
			I0(1) = GF0(1,0)*(r_V(0,0)-r_V(0,1)) + GF0(1,1)*(r_V(1,0)-r_V(1,1)) - G0F0(1,0)*r_V(0,1) - G0F0(1,1)*r_V(1,1);

			I3(0) = GF1(0,0)*(r_V(0,2)-r_V(0,3)) + GF1(0,1)*(r_V(1,2)-r_V(1,3)) - G0F1(0,0)*r_V(0,3) - G0F1(0,1)*r_V(1,3);
			I3(1) = GF1(1,0)*(r_V(0,2)-r_V(0,3)) + GF1(1,1)*(r_V(1,2)-r_V(1,3)) - G0F1(1,0)*r_V(0,3) - G0F1(1,1)*r_V(1,3);

			double SIE = abs(I3(1)/I0(0));
			double Ic_real;
			if (igate_ptr->gatetype != "PRI_INPUT")
			{	
				Ic_real = Vdd/(rho_F*igate_ptr->optM_T/(0.5*igate_ptr->optM_W*igate_ptr->optM_L) + rho_N*ogate_ptr->iptC_T[ipt_order]/(0.5*igate_ptr->optM_W*igate_ptr->optM_L));
			}
			else if (igate_ptr->gatetype == "PRI_INPUT")
			{	
				Ic_real = Vdd/(rho_F*igate_ptr->optM_T/(igate_ptr->optM_W*igate_ptr->optM_L) + rho_N*ogate_ptr->iptC_T[ipt_order]/(igate_ptr->optM_W*igate_ptr->optM_L));
			}
	
			double Ns = Ms*ogate_ptr->optM_L*ogate_ptr->optM_W*ogate_ptr->optM_T/mu_B;
			double t_sw = 2*fsw*q*Ns/abs(Ic_real*SIE);
		//	cout<<t_sw<<endl;
		
			ogate_ptr->ext_ipt_delay[ipt_order] = t_sw;
		//	cout<<ogate_ptr->ext_ipt_delay[ipt_order]<<endl;

		}

	}
	else if ( ( std::find(gatetypeconstlist.begin(),gatetypeconstlist.end(),ogate_ptr->gatetype) != gatetypeconstlist.end() ) || ( std::find(gatetypepriinputlist.begin(),gatetypepriinputlist.end(),ogate_ptr->gatetype) != gatetypepriinputlist.end() ) )
	{
		cout<<"Error: This Output Gate Should Not Enter This Fuction For No Input To It!!!"<<endl;
		assert(false);
	}
	else
	{
		cout<<"Not Defined Output GateType: "<<ogate_ptr->gatetype<<endl;
		assert(false);
	}

}
