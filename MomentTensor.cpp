#include "ModelInfo.h"
#include <iostream>

int main(int argc, char* argv[]) {
	/* check #params */
	if( argc!=5 && argc!=6 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [strike] [dip] [rake] [type(0=NED, 1=USE)] [M0(optional)]"<<std::endl;
		exit(-1);
	}

	// take input
	bool USE = atoi(argv[4])==1;
	float M0 = argc==6 ? atof(argv[5]) : 1.;

	// check focal info
	FocalInfo<float> finfo( atof(argv[1]), atof(argv[2]), atof(argv[3]), 0. );
	if( ! finfo.isValid() ) {
		std::cerr<<"Invalid input focal info: "<<finfo<<std::endl;
		exit(-2);
	}

	// output
	finfo.printMomentTensor(M0, USE);
	finfo.printMTDecomposed(M0, USE);

	return 0;
}
