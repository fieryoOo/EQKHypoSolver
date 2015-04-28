#include "DataTypes.h"
#include "EQKAnalyzer.h"
#include "VectorOperations.h"
#include "SDContainer.h"
#include "ModelSpace.h"
#include "Searcher.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

int main( int argc, char* argv[] ) {
	if( argc != 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [param file]"<<std::endl;
		return -1;
	}

	try {

		// ********** Preparations ********** //
		// initialize model space
		ModelSpace ms( argv[1] );
		// initialize eqk analyzer
		EQKAnalyzer eka( argv[1] );
		eka.LoadData();

		// ********** iterative simulated annealing ********** //
		// search for epicenter and focal mechanism separately
		int niter = 3;
		ms.SetFreeFocal();	// allow perturbing to any focal mechanism, but start at the input focal info
		int nsearch = 16384, Tfactor = 16;
		for( int iter=0; iter<niter; iter++ ) {
			// search for epicenter
			ms.FixFocal();				// have focal mechanism fixed
			eka.PredictAll( ms );	// not necessary, but following search runs faster since Focal is fixed
			if( iter==0 ) eka.SetInitSearch( true );			// use Love group data only!
			Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, 500, 0., 0, std::cout, false );
			if( iter==0 ) eka.SetInitSearch( false );	// use all data
			// search for focal info
			ms.FixEpic();		// have epicenter fixed
			eka.PredictAll( ms );	// not necessary, but following search runs faster since Epic is fixed
			float alpha = std::pow(0.01/Tfactor,1.25/nsearch);	// alpha is emperically decided
			Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, nsearch, alpha, Tfactor, std::cout, false );	// do not save search info
			nsearch /= 2, Tfactor /= 2;
		}
		ms.unFix();			// free both to perturb

		// ********** monte carlo for posterior distributions ********** //
		ms.BoundFocal();	// constrain focal mechanism to perturb near the current Mstate
		// initial MC search around the SA result
		// search for 3000 steps to stablize
		nsearch = 3000;
		auto SIV = Searcher::MonteCarlo<ModelInfo>( ms, eka, nsearch, std::cout );
		// decide perturb step length for each parameter based on the model sensitivity to them
		// perturb steps are defined to be (ub-lb) * Sfactor, where ub and lb are the boundaries decided by:
		// assuming current model state to be the best fitting model, move away
		// from this state until the probability of acceptance <= Pthreshold
		ms.EstimatePerturbs( eka );
		// second (final) Monte Carlo Search with desired perturb sizes
		nsearch = 30000; SIV.clear();
		SIV = Searcher::MonteCarlo<ModelInfo>( ms, eka, nsearch, std::cout );
		std::ofstream fout("MonteCarlo.log");
		for( auto si : SIV )
			fout << si << "\n";

	} catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      return -1;
   } catch(...) {
      std::cerr<<"Unknown Exception!"<<std::endl;
      return -2;
   }

	return 0;
}
