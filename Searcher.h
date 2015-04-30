#ifndef SEARCHER_H
#define SEARCHER_H

#include "MyOMP.h"
#include <iostream>
#include <chrono>
#include <random>
#include <algorithm>
#include <functional>

namespace Searcher {
	// Interfaces. Required by the searcher!!
	template < class MI >
	class IModelSpace {	
		virtual void Perturb( MI& ) const = 0;
		virtual void SetMState( const MI& ) = 0;
		// also: ModelSpace has to be assignable to MI
	};

	template < class MI >
	class IDataHandler {
		virtual float Energy( const MI&, int& Ndata ) const = 0;
	};


	// to save search information
	template < class MI >
	struct SearchInfo {
		int isearch;	// search#
		int ithread;	// thread#
		float T;			// current temperature
		MI info;			// model info
		int Ndata;		// #data
		float E;			// energy of current model
		int accepted;	// accepted in search?

		SearchInfo() {}
		SearchInfo( int isearch, float T, MI info, int Ndata, float E, bool accepted = false )
			: isearch(isearch), T(T), info(info), Ndata(Ndata), E(E), accepted(accepted), ithread(omp_get_thread_num()) {}

		friend bool operator< ( const SearchInfo& s1, const SearchInfo& s2 ) { return s1.isearch<s2.isearch; }
		friend std::ostream& operator<<( std::ostream& o, const SearchInfo<MI>& mi ) {
			o << "search# " << std::setw(3) << mi.isearch << " (T=" << mi.T << ", tid="<<std::setw(2)<<mi.ithread
			  <<"):\tminfo = (" << mi.info << ")\tN = " << mi.Ndata << "\tE = " << mi.E;
			if( mi.accepted == 1 ) o << " (accepted)";
			else if( mi.accepted == 0 ) o << " (rejected)";
			else o << " (best)";
			return o;
		}
	};


   // accept (likelihood function)
   static bool Accept(const float E, const float Enew, const double T, const float probability) {
      bool accept;
      if( Enew<E || probability<exp((E-Enew)/T) ) {
         accept = true;
      } else {
         accept = false;
      }
      //std::cerr<<"E="<<E<<" Enew="<<Enew<<"  accept="<<accept<<" "<<ftmp<<" "<<exp((E-Enew)/T)<<std::endl;
      return accept;
   }


	// Simulated Allealing. Three classes required for 1) modelinfo, 2) modelspace, and 3) datahandler
	// a vector of search info is returned upon request
	// Tinit: initial temperature ( controls initial temperature. called when Tinit is real )
	// Tfactor: T_initial = E_initial * Tfactor ( controls initial temperature. called when Tfactor is int )
	// alpha: T_current = T_last * alpha ( controls rate of temperature decay )
	// saveSI: -1=no, 0=all, 1=accepted
	template < class MI, class MS, class DH >
	std::vector< SearchInfo<MI> > SimulatedAnnealing( MS& ms, DH& dh, const int nsearch,
																	  const float alpha, const int Tfactor,
																	  std::ostream& sout = std::cout, short saveSI = -1 ) {
		return SimulatedAnnealing<MI>( ms, dh, nsearch, alpha, -(float)Tfactor, sout, saveSI );
	}
	template < class MI, class MS, class DH >
	std::vector< SearchInfo<MI> > SimulatedAnnealing( MS& ms, DH& dh, const int nsearch,
																	  const float alpha, const float Tinit,
																	  std::ostream& sout = std::cout, short saveSI = -1 ) {
		// initialize random number generator
		std::default_random_engine generator1( std::chrono::system_clock::now().time_since_epoch().count() + std::random_device{}() );
		std::uniform_real_distribution<float> d_uniform(0., 1.);
		//std::normal_distribution<float> d_normal(0., 1.);
		auto rnd = std::bind( d_uniform, generator1 );

		// force MS, and DH to be derived from the provided interfaces at compile time
		const IModelSpace<MI>& ims = ms;
		const IDataHandler<MI>& idh = dh;

		bool saveacc = saveSI>=0;
		bool saverej = saveSI==0;
		// initial energy (force MS to be assignable to MI at compile time)
		int Ndata;
		float E = dh.Energy( ms, Ndata ), Ebest = E;
		// initial temperature
		float T = Tinit>0. ? Tinit : std::fabs(E*Tinit);
		// SearchInfo vector
		std::vector< SearchInfo<MI> > VSinfo; 
		VSinfo.reserve( nsearch/2*(saveacc+saverej) );
		// save initial
		SearchInfo<MI> sibest( 0, T, ms, Ndata, E, true );
		if( saveacc ) VSinfo.push_back( sibest );
		sout<<sibest<<"\n";
		// main loop
		#pragma omp parallel for schedule(dynamic, 1) private(Ndata)
		for( int i=0; i<nsearch; i++ ) {
			MI minew;
			ms.Perturb( minew );
			float Enew = dh.Energy( minew, Ndata );
			SearchInfo<MI> si( i+1, T, minew, Ndata, Enew );
			#pragma omp critical
			{ // critical begins
			//if( Enew<E || rnd()<exp((E-Enew)/T) ) {
			if( Accept(E, Enew, T, rnd()) ) {
				// update Energy and model state
				ms.SetMState( minew ); E = Enew;
				// mark as accepted
				si.accepted = true;
				// update (if is) the best
				if( Enew < Ebest ) {
					Ebest = Enew;
					sibest = si;
				}
				if( saveacc ) VSinfo.push_back( si );
			} else if( saverej ) VSinfo.push_back( si );
			// temperature decrease
			T *= alpha;
			// output search info
			sout<<si<<"\n";
			} // critical ends
		}
		// output final result
		sibest.accepted = -1;	// -1 for best fitting model
		sout<<sibest<<"\n";
		std::sort( VSinfo.begin(), VSinfo.end() );
		// set model state to the best fitting model
		ms.SetMState( sibest.info );
		return VSinfo;
	}


	// Monte Carlo search: simulated annealing with 1) temperature fixed at 2 and 2) all search info returned
	template < class MI, class MS, class DH >
	std::vector< SearchInfo<MI> > MonteCarlo( MS& ms, DH& dh, const int nsearch, std::ostream& sout = std::cout ) {
		std::cout<<"### Monte Carlo search started. (#search = "<<nsearch<<") ###"<<std::endl;
		return SimulatedAnnealing<MI>( ms, dh, nsearch, 1., 2.0f, sout, 0 );		// save all search info
	}

	template < class MI, class MS, class DH >
	std::vector< SearchInfo<MI> > MonteCarlo( MS& ms, DH& dh, const int nsearch, const std::string& outname ) {
		std::cout<<"### Monte Carlo search in process... (#search="<<nsearch<<") ###\n"
					<<"### results being streamed into file "<<outname<<" ###"<<std::endl;
		std::ofstream fout( outname, std::ofstream::app );
		return SimulatedAnnealing<MI>( ms, dh, nsearch, 1., 2.0f, fout, -1 );	// do not save search info
	}
}

#endif