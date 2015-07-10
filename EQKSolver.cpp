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

/********** options **********
i: Stop after outputing the initial fits, misfits, 
	waveforms(when doing waveform fitting), chi-square(c), and rad-pattern(r)
c: Computes/outputs the initial chi-square
r: Output the initial source radiation pattern
s: do a single big SA for all params instead of the iterative SA
e: Estimate/print perturbation steps at the initial state
d: Debug mode. Computes the initial chi-square of 
	multiple model states taken in from the standard input.
******************************/

inline float Alpha( const int nsearch, const float Tfactor ) {
	return std::pow(0.01/Tfactor,1.25/nsearch);	// emperically decided alpha
}
int main( int argc, char* argv[] ) {
	if( argc < 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [param file] "
					<<"[options (-i=initial-only  -c=chiSquare-init -s=single-SA -r=rad-patterns -e=estimate-perturb -d=debug-mode)]"<<std::endl;
		return -1;
	}

	try {
		// ********** options ********** //
		std::vector<char> options;
		std::string fparam;
		for( int iarg=1; iarg<argc; iarg++ ) {
			std::string str(argv[iarg]);
			if( str[0] == '-' ) {
				for(int istr=1; istr<str.size(); istr++)
					options.push_back( str[istr] );
			} else {
				if( ! fparam.empty() )
					throw std::runtime_error(" ill-formed input! ");
				fparam = str;
			}
		}

		// ********** Preparations ********** //
		// initialize model space
		ModelSpace ms( fparam ); //ms.M0 = 1.3e23;
		// initialize eqk analyzer
		EQKAnalyzer eka( fparam );
		eka.LoadData();

		// initial output
		eka.OutputFits( ms );
		eka.OutputMisfits( ms );
		eka.OutputWaveforms( ms, eka.outdir_sac + "_init" );

		// option -d: debug mode
		if( std::find(options.begin(), options.end(), 'd') != options.end() ) {
			for(int imodel=1; true; imodel++) {
				std::cout<<"input model state "<<imodel<<" (lon, lat, t0, stk, dip, rak, dep, M0; '#' to stop): ";
				std::string line; std::getline(std::cin, line);
				if( line == "#" ) {
					std::cout<<"  end of input\n";
					break;
				}
				try {
					ModelInfo mi(line);
					float chiS; int Ndata;
					//eka.chiSquare( mi, chiS, Ndata );
					eka.Energy( mi, chiS, Ndata );
					std::cout<<"  chi square = "<<chiS<<" Ndata = "<<Ndata<<std::endl;
				} catch(...) {
					std::cout<<"  What was that? Try again please\n";
				}
			}
			return 0;
		}

		// option -c: print out initial chiSquare
		if( std::find(options.begin(), options.end(), 'c') != options.end() ) {
			float chiS; int Ndata;
			eka.chiSquare( ms, chiS, Ndata );
			std::ofstream fout( eka.outname_misAll, std::ofstream::app );
			if(fout) fout<<"### chiS="<<chiS<<" Ndata="<<Ndata<<" at ("<<static_cast<ModelInfo&>(ms)<<") ###\n";
			float E; eka.Energy(ms, E, Ndata);
			std::cout<<"### chiS="<<chiS<<" Ndata="<<Ndata<<"   E="<<E<<" at ("<<static_cast<ModelInfo&>(ms)<<") ###\n";
		}

		// option -r: output initial source patterns
		if( std::find(options.begin(), options.end(), 'r') != options.end() )
			eka.OutputSourcePatterns(ms);

		// option -e: Estimate/print perturbation steps at the initial state
		if( std::find(options.begin(), options.end(), 'e') != options.end() )
			ms.EstimatePerturbs( eka, 0.15 );

		// option -i: stop after output initial fit and misfits
		if( std::find(options.begin(), options.end(), 'i') != options.end() )
			return 0;

		// option -s: do a single big SA for all params instead of the iterative SA
		bool singleSA = std::find(options.begin(), options.end(), 's') != options.end();

		// ********** Initial simulated annealing to approach global optimum ********** //
		ms.SetFreeFocal();	// allow perturbing to any focal mechanism, but start at the input focal info
		if( singleSA ) {
			int nsearch = 9000, Tfactor = 8;
			float alpha = Alpha(nsearch, Tfactor); 
			auto SIV = Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, nsearch, alpha, Tfactor, std::cout, -1 );	// do not save Sinfo
			VO::Output( SIV, eka.outname_misL, true );	// append to file
			//ms.Centralize(); being called later in Bound()
			// output
			eka.OutputFits( ms );
			eka.OutputMisfits( ms );
		} else {
			// ********** iterative simulated annealing ********** //
			// search for epicenter and focal mechanism separately
			//int niterSA = 3, nsearch = 8192, Tfactor = 16;
			int niterSA = 3, nsearch = 4096, Tfactor = 8;
			for( int iter=0; iter<niterSA; iter++ ) {
				// search for epicenter
				ms.FixFocal();				// have focal mechanism fixed
				eka.PredictAll( ms );	// not necessary, but following search runs faster since Focal is fixed
				if( iter==0 ) eka.SetInitSearch( true );			// use Love group data only!
				auto SIV = Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, 500, 0., 0, std::cout, 1 );
				VO::Output( SIV, eka.outname_misL, true );	// append to file
				if( iter==0 ) eka.SetInitSearch( false );	// use all data
				// search for focal info
				ms.FixEpic();		// have epicenter fixed
				eka.PredictAll( ms );	// not necessary, but following search runs faster since Epic is fixed
				float alpha = Alpha(nsearch, Tfactor);
				SIV = Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, nsearch, alpha, Tfactor, std::cout, 1 );	// save info fo accepted searches
				VO::Output( SIV, eka.outname_misF, true );	// append to file
				// centralize the model space around the current MState
				ms.Centralize();
				// output
				eka.OutputFits( ms );
				eka.OutputMisfits( ms );
				nsearch /= 2, Tfactor /= 2;
			}
			//ms.unFix();	// free both to perturb // not necessary, freed in 'Bound()'
		}

		// ********** Post simulated annealing for further optimization ********** //
		// constrain model to perturb near the current Mstate ( Rparam = ? * (0.15, 0.15, 2, 30, 20, 30, 5) )
		// with a small pertfactor to approach the optimum solution faster
		ms.Bound( 2.5, 0.03 );
		// initial MC search around the SA result to stablize
		int nsearch = 5000, Tfactor = 2;
		//auto SIV = Searcher::MonteCarlo<ModelInfo>( ms, eka, nsearch, std::cout );
		//Searcher::MonteCarlo<ModelInfo>( ms, eka, nsearch, eka.outname_pos );
		float alpha = Alpha(nsearch, Tfactor);
		Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, nsearch, alpha, Tfactor, std::cout, -1 );	// do not save Sinfo
		eka.OutputFits( ms );
		eka.OutputMisfits( ms );

		// ********** monte carlo for posterior distributions ********** //
		// constrain model to perturb near the current Mstate ( Rparam = ? * (0.15, 0.15, 2, 30, 20, 30, 5) )
		// perturbation steps are decided later by EstimatePerturbs
		ms.Bound( 2. );	// set Rfactor = 2.0 to be safe
		// decide perturb step length for each parameter based on the model sensitivity to them
		// perturb steps are defined to be (ub-lb) * sfactor, where ub and lb are the boundaries decided by:
		// assuming current model state to be the best fitting model, move away
		// from this state until the probability of acceptance <= Pthreshold
		ms.EstimatePerturbs( eka, 0.15 );	// sfactor default= 0.1
		// second (final) Monte Carlo Search with desired perturb sizes
		nsearch = 50000; 
		Searcher::MonteCarlo<ModelInfo>( ms, eka, nsearch, eka.outname_pos );

		// final output
		eka.OutputFits( ms );
		eka.OutputMisfits( ms );
		eka.OutputSourcePatterns( ms );
		eka.OutputWaveforms( ms );

	} catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      return -1;
   } catch(...) {
      std::cerr<<"Unknown Exception!"<<std::endl;
      return -2;
   }

	return 0;
}
