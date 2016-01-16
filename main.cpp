#include "gate.h"

int main(int argc, char* argv[])
{
	clock_t start, middle1, middle2, end;
	start = clock();

        if (argc<4)
	{
                cout<<"No enough file name given!"<<endl;
		cout<<"./sta 10nm_ASL.genlib CXXX_mapped.blif CXXXoutxxxW.pl EDoutputfile.txt"<<endl;
		return -1;
        } 
	
	ifstream lib(argv[1]);
	if(!lib){
		cout<<"Cannot open lib file!"<<endl;
		cout<<"./sta 10nm_ASL.genlib CXXX_mapped.blif CXXXoutxxxW.pl EDoutputfile.txt"<<endl;
		assert(false);
	}


	ifstream bliff(argv[2]);
	if(!bliff){
		cout<<"Cannot open blif file!"<<endl;
		cout<<"./sta 10nm_ASL.genlib CXXX_mapped.blif CXXXoutxxxW.pl EDoutputfile.txt"<<endl;
		assert(false);
	}

	ifstream plf(argv[3]);
	if(!plf)
	{
		cout<<"Cannot open pl(placement) file!"<<endl;
		cout<<"./sta 10nm_ASL.genlib CXXX_mapped.blif CXXXoutxxxW.pl EDoutputfile.txt"<<endl;
		assert(false);
	}

	ofstream ofile;
	ofile.open(argv[4]);
	if(ofile.is_open() == 0)
	{
		cout<<"Cannot open Energy-Delay output file!"<<endl;
		cout<<"./sta 10nm_ASL.genlib CXXX_mapped.blif CXXXoutxxxW.pl EDoutputfile.txt"<<endl;
		assert(false);
	}	

// starting to read files
	
	Read_Blif_N_Pl(bliff, plf);
	
	Initialization(ofile);
	
	middle1 = clock();
	cout<<"Preprocess Run Time: "<<(double)(middle1-start)/CLOCKS_PER_SEC<<"S"<<endl;

//	TestOptimization(ofile);

//	PriOptimization(ofile);

	Optimization(ofile);

	middle2 = clock();
	cout<<"Optimization Run Time: "<<(double)(middle2-middle1)/CLOCKS_PER_SEC<<"S"<<endl;

	RecalcWithRounding(ofile);

	PostOperation();
	ofile.close();

	end = clock();
	cout<<"Program Run Time: "<<(double)(end-start)/CLOCKS_PER_SEC<<"S"<<endl;

        return 0;

}
