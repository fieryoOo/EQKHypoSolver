#include "ModelInfo.h"
#include <iostream>

int main(int argc, char* argv[]) {
	/* check #params */
	if( argc!=5 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [strike] [dip] [rake] [M0]"<<std::endl;
		exit(-1);
	}

	/* check input focal info */
	FocalInfo<float> finfo( atof(argv[1]), atof(argv[2]), atof(argv[3]), 0. );
	if( ! finfo.isValid() ) {
		std::cerr<<"Invalid input focal info: "<<finfo<<std::endl;
		exit(-2);
	}

	/* transfer to the auxiliary of the current finfo */
	finfo.MomentTensor( atof(argv[4]) );

	return 0;
}
