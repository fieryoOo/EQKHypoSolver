#include "ModelInfo.h"
#include <iostream>

int main(int argc, char* argv[]) {
	/* check #params */
	if( argc!=4 && argc!=5 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [strike] [dip] [rake] [depth(optional)]"<<std::endl;
		exit(-1);
	}

	/* check input focal info */
	float depth = argc==5 ? atof(argv[4]) : 0.;
	FocalInfo<float> finfo( atof(argv[1]), atof(argv[2]), atof(argv[3]), depth );
	if( ! finfo.isValid() ) {
		std::cerr<<"Invalid input focal info: "<<finfo<<std::endl;
		exit(-2);
	}

	/* transfer to the auxiliary of the current finfo */
	finfo.Auxiliary();
	std::cout<<finfo<<std::endl;

	return 0;
}
