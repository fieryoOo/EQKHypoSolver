#include "DataTypes.h"
#include "EQKAnalyzer.h"
#include "VectorOperations.h"
#include "SDContainer.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

int main( int argc, char* argv[] ) {
	if( argc != 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [param file]"<<std::endl;
		return -1;
	}

	try {
		EQKAnalyzer eka( argv[1] );
		eka.LoadData();
		eka.InitRadPattern();
		int N;
		float chiS, wSum;
std::cerr<<"1\n";
		eka.chiSquare( eka._model, chiS, wSum, N );
std::cerr<<"2\n";
	} catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      return -1;
   } catch(...) {
      std::cerr<<"Unknown Exception!"<<std::endl;
      return -2;
   }

	return 0;
}
