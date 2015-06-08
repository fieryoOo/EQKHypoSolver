#include "SynGenerator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

/* FORTRAN entrance */
#ifdef __cplusplus
extern"C" {
#endif
	void surfsyn_( const char* fparam, const char* feigen, const char* fout,
						const char* wavetype, const char* mode, const float* depth);
/*
   void rad_pattern_r_(char *feig_buff, int *eig_len, int *phvnper, int *phvdper,
                       const float *strike, const float *dip, const float *rake, const float *depth,
                       const float *per, int *nper,
                       float *azi, float grT[][nazi], float phT[][nazi], float amp[][nazi]);
*/
#ifdef __cplusplus
}
#endif


void SynGenerator::Run() {
	char fparam[255] = "RUNR", feigen[255] = "R.R", fout[255] = "bred", wtype[2] = "R", mode[2] = "0";
	std::fill(fparam+4, fparam+255, ' ');
	std::fill(feigen+3, feigen+255, ' ');
	std::fill(fout+4, fout+255, ' ');
	wtype[1] = ' '; mode[1] = ' ';
	float depth = 5.6;
	surfsyn_( fparam, feigen, fout, wtype, mode, &depth );
}

bool SynGenerator::Synthetic( const float lon, const float lat, const std::string& chname, 
										const float f1, const float f2, const float f3, const float f4, SacRec& sac ) {
	std::string saclist("sacS.list");
	bool found = false;
	std::ifstream fin(saclist);
	if( ! fin )
		throw std::runtime_error("cannot access file "+saclist);
	for( std::string line; std::getline(fin, line); ) {
		std::stringstream ss(line); ss >> line;
		sac.Load( line );
		if( chname==sac.chname() && lon==sac.shd.stlo && lat==sac.shd.stla ) {
			found = true;
			break;
		}
	}
	if( ! found ) {
		sac.clear();
	} else {
		sac.Mul(-1.);
		sac.Filter(f1,f2,f3,f4);
		sac.shd.b -= minfo.t0;
	}
	return found;
}