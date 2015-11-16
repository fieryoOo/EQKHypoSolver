#include "RadPattern.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

int main( int argc, char* argv[] ) {
	if( argc != 10 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [feigen] [fphvel] [strike] [dip] [rake] [dep] [perlst] [R/L] [outname]"<<std::endl;
		return -1;
	}

	std::string feigname(argv[1]), fphvname(argv[2]);
	float stk = atof(argv[3]), dip = atof(argv[4]), rak = atof(argv[5]), dep = atof(argv[6]);
	char type = argv[8][0];

	try {
		std::vector<float> perlst;
		std::ifstream fin( argv[7] );
		if( ! fin )
			throw std::runtime_error( std::string("cannot read from file ") + argv[7] );
		for( std::string line; std::getline(fin, line); ) {
			float per;
			if( sscanf(line.c_str(), "%f", &per)!=1 || per<=0. ) continue;
			perlst.push_back( per );
		}

		RadPattern rp;
		rp.Predict( type, feigname, fphvname, stk, dip, rak, dep, 1.0, perlst );
		// get source predictions and output
		rp.OutputPreds( argv[9] );
	} catch ( std::exception& e ) {
		std::cerr<<"Exception: "<<e.what()<<std::endl;
	} catch (...) {
		std::cerr<<"Unknown exception!"<<std::endl;
	}

	return 0;
}
