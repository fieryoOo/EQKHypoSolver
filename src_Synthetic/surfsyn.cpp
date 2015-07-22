#include "SynGenerator.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>

int main( int argc, char* argv[] ) {
	// input params
	if( argc != 6 && argc != 7 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [fmodel] [feigen] [wavetype] [mode# (0=fundamental)] [real sac list] [fix_vel (optional)]"<<std::endl;
		return -1;
	}
	// mode #
	int mode = atoi(argv[4]);
	if( mode != atof(argv[4]) )
		throw std::runtime_error("invalid mode# input (expecting integer)");
	/*/ depth
	float depth = atof(argv[5]);
	if( depth<0 || depth>100 )
		throw std::runtime_error("invalid depth input (expecting 0.0-100.0)");
	*/
	if( argc==7 ) {
		float vel = atof(argv[6]);
		if( vel<0.3 || vel>10. ) 
			throw std::runtime_error(std::string("invalid fix_vel input: ")+argv[5]);
	}

	// construct SynGenerator object
	std::string name_fphvel = std::string(argv[2]) + ".phv";
	SynGenerator synG( argv[1], name_fphvel, argv[2], argv[3][0], mode );

	// station list
	//std::string name_fsta("../data/Station.list");
	//synG.LoadSta( name_fsta );
	synG.ClearSta();
	//std::ifstream fin("sac_real.list");
	std::ifstream fin(argv[5]);
	if( ! fin )
		throw std::runtime_error( std::string("IO failed on ")+argv[5] );
	for( std::string line; std::getline(fin, line); ) {
		SacRec sac(line); sac.LoadHD();
		synG.PushbackSta( sac );
	}

	//ModelInfo mi( -114.9049, 41.1453, 0., 38.394, 61.434, -112.444, 6.0, 1.04e23 );
	ModelInfo mi( 245.0974, 41.1448, 0.3254,   251.735, 34.119, -51.516, 6.381, 1.144e23 );

	synG.SetEvent( mi );

//	int npts = 5001;
	int npts = 512;
	float delta = 1.0, perl = 8., perh = 18.;
	float f1 = 0.8/perh, f2 = 1./perh, f3 = 1./perl, f4 = 1.2/perl;
	fin.clear(); fin.seekg(0);
	for( std::string line; std::getline(fin, line); ) {
		SacRec sac(line); sac.LoadHD();
		SacRec sacz, sacn, sace;
		synG.ComputeSyn( sac.stname(), sac.shd.stlo, sac.shd.stla, npts, delta, f1, f2, f3, f4, sacz, sacn, sace, false );
		// write seismograms
		auto wsac = [&]( SacRec& sac ) {
			std::string outname = "sac/"+sac.stname()+"."+sac.chname()+".SAC";
			sac.Write(outname);
		};
		wsac(sacz); wsac(sacn); wsac(sace);
	}

	return 0;
}
