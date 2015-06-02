#include "DataTypes.h"
#include "EQKAnalyzer.h"
#include "VectorOperations.h"
#include "SDContainer.h"
#include "ModelSpace.h"
#include "Searcher.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>

int main( int argc, char* argv[] ) {
	if( argc < 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [param file] "<<std::endl;
		return -1;
	}

	try {
		// ********** options ********** //
		std::string fparam( argv[1] );

		// ********** Preparations ********** //
		// initialize model space
		ModelSpace ms( fparam );
		// initialize eqk analyzer
		EQKAnalyzer eka( fparam );
		eka.LoadData();

		// compute
		float chiS;
		int Ndata;
		eka.Output(ms);
		eka.OutputSourcePatterns(ms);

		eka.chiSquare( ms, chiS, Ndata );
		std::cout<<ms<<"\n"<<chiS<<" "<<Ndata<<std::endl;

	} catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      return -1;
   } catch(...) {
      std::cerr<<"Unknown Exception!"<<std::endl;
      return -2;
   }

}
