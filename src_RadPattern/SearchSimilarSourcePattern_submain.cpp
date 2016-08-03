#include "RadPattern.h"
#include "ModelSpace.h"
#include "Searcher.h"
#include "MyOMP.h"
#include "Timer.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

//struct ModelInfo { float stk, dip, rak, dep; };

class RadPatternDiff : public Searcher::IDataHandler<ModelInfo> {
public:
	//RadPatternDiff() {}

	void InitRayl( const std::string& feigname, const MA3 &sigmasM ) {
		sigmasMR = sigmasM; rp1R.SetModel( 'R', feigname );
		rp2RV.assign(omp_get_max_threads(), rp1R); // to save redundant readings of feig and fphv
		perlstR.clear(); for( auto &paPair : sigmasM ) perlstR.push_back(paPair.first);
	}

	void InitLove( const std::string& feigname, const MA3 &sigmasM ) {
		sigmasML = sigmasM; rp1L.SetModel( 'L', feigname );
		rp2LV.assign(omp_get_max_threads(), rp1L); // to save redundant readings of feig and fphv
		perlstL.clear(); for( auto &paPair : sigmasM ) perlstL.push_back(paPair.first);
	}
	/*
	void SetSigmasL( std::map<float, std::array<float,3>> sigmasM ) {
		sigmasML = sigmasM; perlstL.clear(); for( auto &paPair : sigmasM ) perlstL.push_back(paPair.first);
	}
	*/

	void setRP1( const ftype stk, const ftype dip, const ftype rak, const ftype dep, const ftype M0 = 1. ) {
	   if(! perlstR.empty()) rp1R.Predict( stk, dip, rak, dep, M0, perlstR ); 
	   if(! perlstL.empty()) rp1L.Predict( stk, dip, rak, dep, M0, perlstL ); 
	}

	void setRP2( const ftype stk, const ftype dip, const ftype rak, const ftype dep, const ftype M0 = 1. ) const {
		auto ithread = omp_get_thread_num();
	   if(! perlstR.empty()) rp2RV[ithread].Predict( stk, dip, rak, dep, M0, perlstR ); 
	   if(! perlstL.empty()) rp2LV[ithread].Predict( stk, dip, rak, dep, M0, perlstL ); 
	}

	void OutputDiff( const std::string& outname, const ftype stk, const ftype dip, const ftype rak, const ftype dep , const ftype M0 = 1. ) const {
		setRP2( stk, dip, rak, dep );
		rp2RV[omp_get_thread_num()].OutputPreds("testR_3");
		(rp1R - rp2RV[omp_get_thread_num()]).OutputPreds("testR_1-3");
	}

	void Energy( const ModelInfo& mi, float& E, int& Ndata ) const {
		setRP2( mi.stk, mi.dip, mi.rak, mi.dep );
		auto ithread = omp_get_thread_num(); 

		auto chiSA = perlstR.empty() ? rp2LV[ithread].chiSquare(rp1L, sigmasML) :
						 perlstL.empty() ? rp2RV[ithread].chiSquare(rp1R, sigmasMR) :
						 chiSquare( rp2RV[ithread], rp2LV[ithread], rp1R, rp1L, sigmasMR, sigmasML );
		E = chiSA[0] + chiSA[1] + chiSA[2]; Ndata = chiSA[3]; mi.M0 = perlstR.empty() ? rp2LV[ithread].M0 : rp2RV[ithread].M0;
		//std::cerr<<E<<" "<<Ndata<<" "<<mi.M0<<std::endl; exit(-3);
	}

private:
	RadPattern rp1R, rp1L;
	mutable std::vector<RadPattern> rp2RV, rp2LV;
	bool useR = false, useL = false;
	MA3 sigmasMR, sigmasML;
	std::vector<float> perlstR, perlstL;
};


int main( int argc, char* argv[] ) {
   if( argc != 11 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [R/L/B] [eigen_file_R (.R)] [eigen_file_L (.L)]";
      std::cerr<<" [per-sigmas_R] [per-sigmas_L] [strike] [dip] [rake] [depth] [out_name]"<<std::endl;
      exit(-1);
   }

   // read in type 
   char type = argv[1][0]; bool useR = true, useL = true;
	if( type == 'R' ) { useL = false; } 
	else if( type == 'L' ) { useR = false; } 
	else if( type == 'B' ) {}
	else {
      std::cerr<<"Unknown type: "<<type<<std::endl;
      exit(0);
   }

	// feigs
	std::string feignmR( argv[2] ), feignmL( argv[3] );

   // read in the per-sigmas lists
	// per sigmaG(sec) sigmaP(sec) sigmaA(fraction 0-1)
	auto loadfSigmas = []( const std::string &fsigmas, MA3 &sigmasM ) {
	   std::ifstream fin(fsigmas);
		if( ! fin ) {
	      std::cerr<<"Error(main): Cannot read from file "<<fsigmas<<std::endl;
		   exit(0);
	   }
	   for(std::string line; std::getline(fin, line); ) {
		   float per, sG, sP, sA;
			if( sscanf(line.c_str(), "%f %f %f %f", &per, &sG, &sP, &sA) != 4 ) continue;
	      sigmasM[per] = {sG, sP, sA};
		}
	   fin.close();
		std::cout<<"### "<<sigmasM.size()<<" periods read in from "<<fsigmas<<" ###"<<std::endl;
	};
	MA3 sigmasMR, sigmasML;
	loadfSigmas( argv[4], sigmasMR ); loadfSigmas( argv[5], sigmasML );

   // read in focal info 
   const float stk0 = atof(argv[6]), dip0 = atof(argv[7]), rak0 = atof(argv[8]), dep0 = atof(argv[9]), M0 = 1.;
   std::cout<<"### Input Focal info = ("<<stk0<<" "<<dip0<<" "<<rak0<<" "<<dep0<<" "<<M0<<"). ###"<<std::endl;

	// initialize RadPatternDiff
	RadPatternDiff rpd;
	if(useR) rpd.InitRayl( feignmR, sigmasMR );
	if(useL) rpd.InitLove( feignmL, sigmasML );
	rpd.setRP1( stk0, dip0, rak0, dep0 );

	// initialize model space, fix lon, lat, t0, and M0
	ModelSpace ms( ModelInfo(0., 0., 0., 270, 35, -16, 4.3) );
	ms.SetFreeFocal(); ms.SetPerturb(false, false, false, false, true, true, true, true);
	//ModelSpace ms2( ModelInfo(0., 0., 0., 275, 39, 0, 4.4) );

	// debug
	/*
	float E; int N; ModelInfo mi;
	mi = {0., 0., 0., 247.8, 33.3, -59.6, 5.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 248.8, 33.3, -59.6, 5.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 249.8, 33.3, -59.6, 5.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 250.8, 33.3, -59.6, 5.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 247.8, 34.3, -59.6, 5.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 247.8, 35.3, -59.6, 5.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 247.8, 36.3, -59.6, 5.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 247.8, 33.3, -60.6, 5.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 247.8, 33.3, -61.6, 5.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 247.8, 33.3, -62.6, 5.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 247.8, 33.3, -59.6, 6.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 247.8, 33.3, -59.6, 7.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	mi = {0., 0., 0., 247.8, 33.3, -59.6, 8.2}; rpd.Energy(mi, E, N); std::cerr<<E<<" "<<N<<" ("<<mi<<")"<<std::endl;
	exit(-3); 
	*/

	//Timer timer;
	float Emean, Estd; Searcher::EStatistic<ModelInfo>(ms, rpd, 1000, Emean, Estd);
	//std::cerr<<Emean<<" "<<Estd<<"   CPUtimer = "<<timer.CPUSec()<<" Walltime = "<<timer.WallSec()<<std::endl;

	// run simulated annealing
	std::ofstream fout(argv[10]);	int nsearch = 10000, niter = 1000;
	for(int i=0; i<niter; i++) {
		ms.RandomState();
		Searcher::SimulatedAnnealing<ModelInfo>(ms, rpd, nsearch, Emean+Estd*3, 0.1, 0, fout);
	}

   return 0;
}

