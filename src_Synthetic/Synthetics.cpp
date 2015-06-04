#include "Synthetics.h"
#include <iostream>
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


void Synthetics::Run() {
	char fparam[255] = "RUNR", feigen[255] = "R.R", fout[255] = "bred", wtype[2] = "R", mode[2] = "0";
	std::fill(fparam+4, fparam+255, ' ');
	std::fill(feigen+3, feigen+255, ' ');
	std::fill(fout+4, fout+255, ' ');
	wtype[1] = ' '; mode[1] = ' ';
	float depth = 5.6;
	surfsyn_( fparam, feigen, fout, wtype, mode, &depth );
}
