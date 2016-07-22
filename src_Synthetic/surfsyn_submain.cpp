#include "SynGenerator.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>

int main( int argc, char* argv[] ) {
	// input params
	if( argc != 7 && argc != 8 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [fmodel] [feigen (the phv file should be named ${feigen}.phv)] [wavetype] [mode# (0=fundamental)] [real sac list] [outdir (has to exist)] [f_ModelInfo]"<<std::endl;
		return -1;
	}
	// mode #
	int mode = atoi(argv[4]);
	if( mode != atof(argv[4]) )
		throw std::runtime_error("invalid mode# input (expecting integer)");
	/*
	if( argc==8 ) {
		float vel = atof(argv[7]);
		if( vel<0.3 || vel>10. ) 
			throw std::runtime_error(std::string("invalid fix_vel input: ")+argv[5]);
	}
	*/

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
	ModelInfo mi;
	if( argc == 8 ) {
		std::ifstream finMI( argv[7] );
		std::string MIline;
		std::getline(finMI, MIline);
		mi = ModelInfo(MIline);
	} else {
		mi = ModelInfo( 245.096, 41.154, 0.0,   248.0, 33.5, -54.0, 5.9, 1.1e23 );
	}

	synG.SetEvent( mi );

//	int npts = 5001;
	int npts = 1024;
	float delta = 1.0, perl = 5., perh = 47.5;
	//float f1 = 0.8/perh, f2 = 1./perh, f3 = 1./perl, f4 = 1.2/perl;
	float f1 = 1./60., f2 = 1./55., f3 = 1./4.5, f4 = 1./4.;
	fin.clear(); fin.seekg(0);
	for( std::string line; std::getline(fin, line); ) {
		std::cerr<<line<<std::endl;
		SacRec sac(line); sac.LoadHD();
		SacRec sacz, sacn, sace;
		//if( ! synG.ComputeSyn( sac.stname(), sac.shd.stlo, sac.shd.stla, npts, delta, sacz, sacn, sace, false, f1, f2, f3, f4 ) )
		if( ! synG.ComputeSyn( sac.stname(), sac.shd.stlo, sac.shd.stla, npts, delta, sacz, sacn, sace, false ) )
			continue;
		// write seismograms
		auto wsac = [&]( SacRec& sac ) {
			std::string outname = std::string(argv[6])+"/"+sac.stname()+"."+sac.chname()+".SAC";
			sac.Write(outname);
		};
		if( synG.type == 'L' ) wsac(sace);
		else { wsac(sacz); wsac(sacn); }
	}

	return 0;
}
