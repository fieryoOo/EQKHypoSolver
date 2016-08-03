#include "RadPattern.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

int main( int argc, char* argv[] ) {
	if( argc != 9 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [feigen] [strike] [dip] [rake] [dep] [perlst] [R/L] [outname]"<<std::endl;
		return -1;
	}

	std::string feigname(argv[1]);
	float stk = atof(argv[2]), dip = atof(argv[3]), rak = atof(argv[4]), dep = atof(argv[5]);
	char type = argv[7][0];

	try {
		std::vector<float> perlst;
		std::ifstream fin( argv[6] );
		if( ! fin )
			throw std::runtime_error( std::string("cannot read from file ") + argv[6] );
		for( std::string line; std::getline(fin, line); ) {
			float per;
			if( sscanf(line.c_str(), "%f", &per)!=1 || per<=0. ) continue;
			perlst.push_back( per );
		}

		RadPattern rp( type, feigname );
		rp.Predict( stk, dip, rak, dep, 1.0, perlst );
		// get source predictions and output
		rp.OutputPreds( argv[8] );
	} catch ( std::exception& e ) {
		std::cerr<<"Exception: "<<e.what()<<std::endl;
	} catch (...) {
		std::cerr<<"Unknown exception!"<<std::endl;
	}

	return 0;
}
