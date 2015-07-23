#include "SacRec.h"
#include "Curve.h"
#include <iostream>
#include <fstream>


int main( int argc, char* argv[] ) {
	if( argc!=3 && argc!=4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [sac1] [sac2] [fout(optional)]"<<std::endl;
		return -1;
	}

	try {
		// load
		SacRec sac1, sac2;// saccor;
		sac1.Load(argv[1]);
		sac2.Load(argv[2]);
		//sac1.CrossCorrelate(sac2, saccor);

		// FFT
		SacRec sac_am, sac_ph1, sac_ph2;
		sac1.ToAmPh(sac_am, sac_ph1);
		sac2.ToAmPh(sac_am, sac_ph2);

		// difference
		sac_ph1.Subf( sac_ph2 );

		// correct 2pi and convert from radian to sec
		const float TWO_PI = M_PI * 2.;
		const float f0 = sac_ph1.shd.b, df = sac_ph1.shd.delta;
		auto convert = [&]( const float f, float& val) {
			if( val >= M_PI ) val -= TWO_PI;
			else if( val < -M_PI ) val += TWO_PI;
			val *= 1.0/(f*TWO_PI);
		};
		sac_ph1.Transform2( convert, 1 );
		sac_ph1.sig[0] = 0.;

		// output
		if( argc == 4 ) sac_ph1.Write(argv[3]);
		else sac_ph1.Dump();

	} catch(...) {
		std::cerr<<"Unknown exception!"<<std::endl;
		exit(-2);
	}

	return 0;
}
