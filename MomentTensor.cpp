#include "ModelInfo.h"
#include <iostream>

int main(int argc, char* argv[]) {
	/* check #params */
	if( argc!=4 && argc!=5 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [strike] [dip] [rake] [M0(optional)]"<<std::endl;
		exit(-1);
	}

	/* check input focal info */
	float M0 = argc==5 ? atof(argv[4]) : 1.;
	FocalInfo<float> finfo( atof(argv[1]), atof(argv[2]), atof(argv[3]), 0. );
	if( ! finfo.isValid() ) {
		std::cerr<<"Invalid input focal info: "<<finfo<<std::endl;
		exit(-2);
	}

	/* transfer to the auxiliary of the current finfo */
	finfo.MomentTensor(M0);

	return 0;
}
