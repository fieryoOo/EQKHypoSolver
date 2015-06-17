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
		std::cerr<<"Usage: "<<argv[0]<<" [param file] "
					<<"[options (-i=initial_only  -c=chiSquare_init -s=skip_SA -r=rad_patterns)]"<<std::endl;
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
		ModelSpace ms( fparam );
		// initialize eqk analyzer
		EQKAnalyzer eka( fparam );
		eka.LoadData();

/*
float chiS; int Ndata;
eka.chiSquareW( ms, chiS, Ndata );
std::cerr<<"chi square = "<<chiS<<" "<<Ndata<<std::endl;
exit(-3);
*/

		// initial output
		eka.Output( ms );
		eka.OutputMisfits( ms );

		// option -c: print out initial chiSquare
		if( std::find(options.begin(), options.end(), 'c') != options.end() ) {
			float chiS; int Ndata;
			eka.chiSquareM( ms, chiS, Ndata );
			std::ofstream fout( eka.outname_misAll, std::ofstream::app );
			if(fout) fout<<"### chiS="<<chiS<<" Ndata="<<Ndata<<" at ("<<static_cast<ModelInfo&>(ms)<<") ###\n";
			float E; eka.Energy(ms, E, Ndata);
			std::cout<<"### chiS="<<chiS<<" Ndata="<<Ndata<<"   E="<<E<<" at ("<<static_cast<ModelInfo&>(ms)<<") ###\n";
		}

		// option -r: output initial source patterns
		if( std::find(options.begin(), options.end(), 'r') != options.end() )
			eka.OutputSourcePatterns(ms);

		// option -i: stop after output initial fit and misfits
		if( std::find(options.begin(), options.end(), 'i') != options.end() )
			return 0;

		// ********** iterative simulated annealing ********** //
		// search for epicenter and focal mechanism separately
		int niter = 3;
		// option -s: skip simulated annealing
		if( std::find(options.begin(), options.end(), 's') != options.end() ) niter = 0;
		ms.SetFreeFocal();	// allow perturbing to any focal mechanism, but start at the input focal info
		int nsearch = 8192, Tfactor = 16;
		for( int iter=0; iter<niter; iter++ ) {
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
			float alpha = std::pow(0.01/Tfactor,1.25/nsearch);	// alpha is emperically decided
			SIV = Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, nsearch, alpha, Tfactor, std::cout, 1 );	// save info fo accepted searches
			VO::Output( SIV, eka.outname_misF, true );	// append to file
			// centralize the model space around the current MState
			ms.Centralize();
			// output
			eka.Output( ms );
			eka.OutputMisfits( ms );
			nsearch /= 2, Tfactor /= 2;
		}
		//ms.unFix();	// free both to perturb // not necessary, freed in 'Bound()'

		// ********** monte carlo for posterior distributions ********** //
		// constrain model to perturb near the current Mstate ( Rparam = ? * (0.15, 0.15, 2, 30, 20, 30, 5) )
		// with a small pertfactor to approach the optimum solution faster
		ms.Bound( 2.5, 0.03 );
		// initial MC search around the SA result to stablize
		nsearch = 5000;
		//auto SIV = Searcher::MonteCarlo<ModelInfo>( ms, eka, nsearch, std::cout );
		//Searcher::MonteCarlo<ModelInfo>( ms, eka, nsearch, eka.outname_pos );
		Tfactor = 2; float alpha = std::pow(0.01/Tfactor,1.25/nsearch);	// alpha is emperically decided
		Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, nsearch, alpha, Tfactor, std::cout, -1 );	// do not save Sinfo
		eka.Output( ms );
		eka.OutputMisfits( ms );
		// decide perturb step length for each parameter based on the model sensitivity to them
		// perturb steps are defined to be (ub-lb) * sfactor, where ub and lb are the boundaries decided by:
		// assuming current model state to be the best fitting model, move away
		// from this state until the probability of acceptance <= Pthreshold
		ms.Bound( 2. );	// set Rfactor = 2.0 to be safe
		ms.EstimatePerturbs( eka, 0.15 );	// sfactor default= 0.1
		// second (final) Monte Carlo Search with desired perturb sizes
		nsearch = 50000; 
		Searcher::MonteCarlo<ModelInfo>( ms, eka, nsearch, eka.outname_pos );

		// final output
		eka.Output( ms );
		eka.OutputMisfits( ms );
		eka.OutputSourcePatterns( ms );

	} catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      return -1;
   } catch(...) {
      std::cerr<<"Unknown Exception!"<<std::endl;
      return -2;
   }

	return 0;
}
