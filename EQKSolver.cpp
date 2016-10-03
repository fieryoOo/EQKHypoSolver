#include "DataTypes.h"
#include "EQKAnalyzer.h"
#include "VectorOperations.h"
#include "SDContainer.h"
#include "ModelSpace.h"
#include "Searcher.h"
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>

/********** options **********
i: Stop after outputing the initial fits, misfits, 
	waveforms(when doing waveform fitting), chi-square(c), and rad-pattern(r)
c: Computes/outputs the initial chi-square
r: Output the initial source radiation pattern
s: do a single big SA for all params instead of the iterative SA
g: run only a fast SA to stablize, followed by the Monte-Carlo search (assuming close-enough input model info)
m: start the Monte-Carlo search immediately (assuming highly-optimized input model info)
n: do not run the Monte-Carlo search
e: Estimate/print perturbation steps at the initial state
d: Debug mode. Computes the initial chi-square of 
	multiple model taken in from the standard input.
******************************/

/* moved to inside Searcher
inline float Alpha( const int nsearch, const float Tfactor ) {
	return std::pow(0.01/Tfactor,1.25/nsearch);	// emperically decided alpha
}
*/
int main( int argc, char* argv[] ) {
	if( argc < 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [param file] "
					<<"[options (-ic=(stop after)initial-computation -io=(stop after)initial-outputs -pic=print-init-chiSquare -dm=debug-mode -oir=output-init-radpatterns "
					<<"-uev=use-estimated-variances -pt?=-perturb-type-(0/1/2) -rsa?-?=regular-SA -lc=linear-cooling -gi=good-initial -mco?=Monte-Carlo-only -nmc=No-Monte_Carlo)]"<<std::endl;
		return -1;
	}
		// option -v: output initial variances across the whole region

	try {
		// ********** options ********** //
		//std::vector<char> options;
		std::unordered_map<std::string, std::string> options;
		std::string fparam;
		for( int iarg=1; iarg<argc; iarg++ ) {
			std::string str(argv[iarg]);
			if( str[0] == '-' ) {
				int inum = str.find_first_of("0123456789");
				std::string arg;
				if(inum == std::string::npos) inum = str.size();
				else arg = str.substr(inum,str.size()-inum);
				options.emplace(str.substr(1,inum-1), arg);
				//for(int istr=1; istr<str.size(); istr++) options.push_back( str[istr] );
			} else {
				if( ! fparam.empty() )
					throw std::runtime_error(" ill-formed input! ");
				fparam = str;
			}
		}

		// ********** Preparations ********** //
		// initialize model space
		ModelSpace ms( fparam ); //ms.M0 = 1.3e23;
		// initialize eqk analyzer (do not move/save old outputs just yet)
		EQKAnalyzer eka( fparam, false );
		eka.LoadData();

		// option -pic: print out initial chiSquare
		//if( std::find(options.begin(), options.end(), 'c') != options.end() ) {
		if( options.find("pic") != options.end() ) {
			float chiS; int Ndata;
			eka.chiSquare( ms, chiS, Ndata );
			std::ofstream fout( eka.outname_misAll, std::ofstream::app );
			//if(fout) fout<<"### chiS="<<chiS<<" Ndata="<<Ndata<<" at ("<<static_cast<ModelInfo&>(ms)<<") ###\n";
			//float E; eka.Energy(ms, E, Ndata);
			std::cout<<"### chiS="<<chiS<<" Ndata="<<Ndata<<"   E="<<chiS*eka._indep_factor<<" at ("<<static_cast<ModelInfo&>(ms)<<") ###\n";
		}

		// option -dm: debug mode
		//if( std::find(options.begin(), options.end(), 'd') != options.end() ) {
		if( options.find("dm") != options.end() ) {
			for(int imodel=1; true; ) {
				std::cout<<"input model "<<imodel<<" (lon, lat, t0, stk, dip, rak, dep, M0; '#' to stop): ";
				std::string line; std::getline(std::cin, line);
				if( line == "#" ) {
					std::cout<<"  end of input\n";
					break;
				}
				try {
					ModelInfo mi(line); imodel++;
					float chiS; int Ndata;
					eka.chiSquare( mi, chiS, Ndata );
					//eka.Energy( mi, E, Ndata );
					std::cout<<"  chi square = "<<chiS<<" Ndata = "<<Ndata<<std::endl;
				} catch(...) {
					std::cout<<"  What was that? Try again please\n";
				}
			}
			return 0;
		}

		// option -ic(initial computation): stop before output initial fit and misfits
		if( options.find("ic") != options.end() ) return 0;

		// initial output
		eka.SaveOldOutputs();
		eka.OutputFits( ms );
		eka.OutputMisfits( ms );
		eka.OutputWaveforms( ms, eka.outdir_sac + "_init" );
		eka.OutputSigmas();
		// option -oir: output initial rad patterns
		if( options.find("oir") != options.end() ) eka.OutputSourcePatterns(ms);

		// option -io(initial output): stop after output initial fit and misfits
		if( options.find("io") != options.end() ) return 0;

		// option -pt (perturbation type):
		// 0 = use init pert sizes throughout, 1 = estimate (after SA) pert sizes for MC, 2 = estimate (before and after SA) pert sizes for both
		int pt = options.find("pt")==options.end() ? 1 : std::stoi(options["pt"]);
		if( pt == 2 ) ms.EstimatePerturbs( eka, 0.15 );

		// option -rsa: do a single big SA for all params instead of the iterative SA
		bool regularSA = options.find("rsa") != options.end();

		// option -gi: run only a fast SA to stablize, followed by the Monte-Carlo search (assuming close-enough input model info)
		// option -mco: start the Monte-Carlo search immediately (assuming a highly-optimized input model state)
		bool doSA1 = options.find("gi")==options.end() && options.find("mco") == options.end();
		bool doSA2 = options.find("mco")==options.end() && !regularSA;

		// option -nmc: do not run the Monte-Carlo search
		bool doMC = options.find("nmc") == options.end();

		// option -lc: use linear cooling schedule for SA
		int cooltype = options.find("lc") != options.end();	// 0 = exponential cooling, 1 = linear cooling

		// estimate Energy statistics
		eka.SetCorrectM0(true);	// M0 will be corrected to produce least chiS
		float Emean, Estd; int nthd = omp_get_max_threads();
		ms.SetFreeFocal();	// allow perturbing to any focal mechanism, but start at the input focal info
		Searcher::EStatistic<ModelInfo>(ms, eka, nthd>10?nthd*5:50, Emean, Estd);
		std::cout<<"### Energy Statistics: Emean = "<<Emean<<"  Estd = "<<Estd<<" ###"<<std::endl;

		// ********** Initialize simulated annealing to approach global optimum ********** //
		double Tinit = (cooltype==0?3:0.5) * (Emean+3*Estd), Tfinal = cooltype==0 ? 1e-4*Emean : 0;
		if( doSA1 ) {

			if( regularSA ) {
				// ********** single simulated annealing ********** //
				// search for epicenter and focal mechanism simultaneously
				//int nsearch = 7200, Tfactor = 8;
				int nsearch = 15000, niter = 1;//, Tfactor = 1000;
				std::stringstream ss(options["rsa"]); std::string tok; 
				if( std::getline(ss,tok,'-') ) nsearch = std::stoi(tok);
				if( std::getline(ss,tok,'-') ) niter = std::stoi(tok);
				//double alpha = Searcher::Alpha(nsearch, Tfactor);
				auto msbest = ms; float Ebest; int Ndata; eka.Energy( ms, Ebest, Ndata );
				for(int iter=0; iter<niter; iter++) {
					ms.SetPerturb( true, true, true, false, true, true, true, true );	// do not perturb M0
					//auto SIV = Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, nsearch, alpha, Tfactor, std::cout, 0, true );	// save info for accepted searches
					auto SIV = Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, nsearch*4/5, Tinit, Tfinal, cooltype, std::cout, 0, true );	// save info for accepted searches
					VO::Output( SIV, eka.outname_misL, true );	// append to file
					if( pt > 0 ) ms.Bound( 2.5, 0.03 );
					SIV = Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, nsearch/5, Tinit*0.01, Tfinal*0.01, cooltype, std::cout, 0, true );	// save info for accepted searches
					VO::Output( SIV, eka.outname_misL, true );	// append to file
					// output
					eka.OutputFits( ms );
					eka.OutputMisfits( ms );
					// update msbest
					if( SIV.back().E < Ebest ) { msbest = ms; Ebest = SIV.back().E; }
					ms.SetFreeFocal(); ms.RandomState();
				}
				ms = msbest;
			} else {
				// ********** iterative simulated annealing ********** //
				// search for epicenter and focal mechanism separately
				int niterSA = 3, nsearch = 8192; //, Tfactor = 16;
				//int niterSA = 2, nsearch = 4096;
				auto Tinit2 = Tinit;
				for( int iter=0; iter<niterSA; iter++ ) {
					// search for epicenter
					ms.SetPerturb( true, true, true, false, false, false, false, false );	// have focal mechanism fixed
					eka.UpdatePreds( ms );	// not necessary, but following search runs faster since Focal is fixed
					if( iter==0 ) eka.SetInitSearch( true );			// use Love group data only!
					auto SIV = Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, 500, 0., 0., 1, std::cout, 0, true );
					VO::Output( SIV, eka.outname_misL, true );	// append to file
					if( iter==0 ) eka.SetInitSearch( false );	// use all data
					// search for focal info
					ms.SetPerturb( false, false, false, false, true, true, true, true );	// have epicenter (and M0) fixed
					eka.UpdatePreds( ms );	// not necessary, but following search runs faster since Epic is fixed
					//double alpha = Searcher::Alpha(nsearch, Tfactor);
					SIV = Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, nsearch, Tinit2, Tfinal, cooltype, std::cout, 0, true );	// save info for accepted searches
					VO::Output( SIV, eka.outname_misF, true );	// append to file
					// centralize the model space around the current MState
					ms.Centralize();
					// output
					eka.OutputFits( ms );
					eka.OutputMisfits( ms );
					nsearch /= 2, Tinit2 /= 2;
				}
				//ms.unFix();	// free both to perturb // not necessary, freed in 'Bound()'
			}
		}

		// ********** secondary simulated annealing for deeper optimization ********** //
		if( doSA2 ) {
			//double Tinit = (cooltype==0?0.01:0.002) * Emean, Tfinal = cooltype==0 ? 1e-7*Emean : 0;
			// constrain model to perturb near the current Mstate ( Rparam = ? * (0.15, 0.15, 2, 30, 20, 30, 5) )
			// with a small pertfactor to approach the optimum solution faster
			if( pt > 0 ) ms.Bound( 2.5, 0.03 );
			ms.SetPerturb( true, true, true, false, true, true, true, true );	// have M0 fixed
			// initial MC search around the SA result to stablize
			int nsearch = 3000;
			//auto SIV = Searcher::MonteCarlo<ModelInfo>( ms, eka, nsearch, std::cout );
			//Searcher::MonteCarlo<ModelInfo>( ms, eka, nsearch, eka.outname_pos );
			//double alpha = Searcher::Alpha(nsearch, Tfactor);
			//ms.SetPerturb( false, false, false, true, false, false, false, false );
			Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, nsearch, Tinit*0.001, Tfinal*0.001, cooltype, std::cout, 0, false );	// do not save Sinfo
			//Searcher::SimulatedAnnealing<ModelInfo>( ms, eka, 10000, alpha, 0.5f, std::cout, -1 );	// do not save Sinfo
			eka.OutputFits( ms );
			eka.OutputMisfits( ms );
			eka.OutputSourcePatterns( ms );
			eka.OutputWaveforms( ms );

		}

		// -uev: estimate/use sigmaS to be the current (after SA) variance across all stations
		if( options.find("uev") != options.end() ) {
			eka.EstimateSigmas( ms ); eka.OutputSigmas();
		}

		if( doMC ) {
			// ********** monte carlo for posterior distributions ********** //
			// constrain model to perturb near the current Mstate ( Rparam = ? * (0.15, 0.15, 2, 30, 20, 30, 5) )
			// perturbation steps are decided later by EstimatePerturbs
			//ms.Bound( 2. );	// set Rfactor = 2.0 to be safe
			eka.SetCorrectM0(false);	// M0 will not be corrected to produce least chiS
			ms.SetFreeFocal();			// allow perturbing to any focal mechanism, but start at the current focal info
			// decide perturb step length for each parameter based on the model sensitivity to them
			// perturb steps are defined to be (ub-lb) * sfactor, where ub and lb are the boundaries decided by:
			// assuming current model to be the best fitting model, move away
			// from this state until the probability of acceptance <= Pthreshold
			if( pt > 0 ) ms.EstimatePerturbs( eka, 0.12 );	// sfactor default= 0.1
			// second (final) Monte Carlo Search with desired perturb sizes
			int nsearch = 100000; 
			if( options.find("mco")!=options.end() && !options["mco"].empty() ) 
				nsearch = std::stoi(options["mco"]);
			Searcher::MonteCarlo<ModelInfo>( ms, eka, nsearch, eka.outname_pos );

			// final output
			eka.OutputFits( ms );				// appended
			eka.OutputMisfits( ms );			// appended
			eka.OutputSourcePatterns( ms );	// overwritten
			eka.OutputWaveforms( ms );			// overwritten
		}

	} catch(std::exception& e) {
      std::cerr<<e.what()<<std::endl;
      return -1;
   } catch(...) {
      std::cerr<<"Unknown Exception!"<<std::endl;
      return -2;
   }

	return 0;
}
