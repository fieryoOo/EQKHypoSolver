#ifndef SEARCHER_H
#define SEARCHER_H

#include "MyOMP.h"
#include "Rand.h"
#include "VectorOperations.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <thread>
#include <random>
#include <algorithm>
#include <functional>

namespace Searcher {
//class Searcher {
	static float _poc = -1.;	// percentage-of-completion of the current search

	// Interfaces. Required by the searcher!!
	template < class MI >
	class IModelSpace {	
		virtual void SetMState( const MI& ) = 0;					// required by SimulatedAnnealing
		virtual void Perturb( MI&, bool isfree ) const = 0;	// and MonteCarlo
		virtual void RandomState( IModelSpace &msnew ) const {}			// required by EStatistic
		// also: IModelSpace has to be assignable to MI
	};

	template < class MI >
	class IDataHandler {
		virtual void Energy( const MI&, float& E, int& Ndata ) const = 0;
	};

	// empirical formula for alpha
	static inline double Alpha( const int nsearch, const double Tfactor ) {
		if( Tfactor <= 0 ) return 0.;
		return std::pow(0.01/Tfactor,1.25/nsearch);	// emperically decided alpha
	}

	// to save search information
	template < class MI >
	struct SearchInfo {
		int isearch;	// search#
		int ithread;	// thread#
		double T;		// current temperature
		MI info;			// model info
		int Ndata = 0;	// #data
		float E = -1.;	// energy of current model
		int accepted;	// accepted in search?

		SearchInfo() {}
		SearchInfo( int isearch, double T, MI info, int Ndata, float E, int accepted = 0 )
			: isearch(isearch), T(T), info(info), Ndata(Ndata), E(E), accepted(accepted), ithread(omp_get_thread_num()) {}

		friend bool operator< ( const SearchInfo& s1, const SearchInfo& s2 ) { return s1.isearch<s2.isearch; }
		friend bool operator== ( const SearchInfo& s1, const SearchInfo& s2 ) { return s1.info==s2.info; }

		friend std::ostream &operator<<( std::ostream& o, const SearchInfo<MI> &si ) {
			o << "search# " << std::setw(3) << si.isearch << " (T=" << si.T << ", tid="<<std::setw(2)<<si.ithread
			  <<"):\tminfo = (" << si.info << ")\tN = " << si.Ndata << "\tE = " << si.E;
			if( si.accepted == 1 ) o << " (accepted)";
			else if( si.accepted == 0 ) o << " (rejected)";
			else if( si.accepted == -1 ) o << " (skipped)";
			else o << " (best)";
			return o;
		}

		friend std::istream &operator>>( std::istream &i, SearchInfo<MI> &si ) {
			std::string stmp;	i >> stmp;
			if( stmp == "search#" ) {
				i>>si.isearch; i.ignore(4,'='); i>>si.T; i.ignore(6,'='); i>>si.ithread; i.ignore(12,'('); 
				i>>si.info; i.ignore(5,'='); i>>si.Ndata; i.ignore(3,'='); i>>si.E; i>>stmp;
				if( stmp == "(accepted)") si.accepted = 1;
				else if( stmp=="(rejected)" ) si.accepted = 0;
				else if( stmp=="(skipped)" ) si.accepted = -1;
				else si.accepted = 2;
			} else {
				i.setstate(std::ios::failbit);
			}
			return i;
		}
	};


   // accept (likelihood function)
   static float Paccept(const float E, const float Enew, const double T) {
		return Enew<=E ? 1. : (T==0?0:exp((E-Enew)/T));
   }


	// estimate Energy n times and compute the statistic
	template < class MI, class MS, class DH >
	void EStatistic( MS& ms, DH& dh, const int n, float &Emean, float &Estd ) {
		std::vector<float> EV(n);
		#pragma omp parallel for schedule(dynamic, 1)
		for( int i=0; i<n; i++ ) {
			MI minew; ms.RandomState( minew );
			int Ndata; dh.Energy( minew, EV[i], Ndata );
		}
		VO::MeanSTD(EV.begin(), EV.end(), Emean, Estd);
		//std::cerr<<Emean<<" "<<Estd<<std::endl;
	} 
	

	// Simulated Annealing. Three classes required for 1) modelinfo, 2) modelspace, and 3) datahandler
	// a vector of search info is returned upon request
	// Tinit: initial temperature ( controls initial temperature. called when Tinit is real )
	// Tfactor: T_initial = E_initial * Tfactor ( controls initial temperature. called when Tfactor is int )
	// alpha: T_current = T_last * alpha ( controls rate of temperature decay )
	// outSI: -1=no, 0=all, 1=accepted
	/*
	template < class MI, class MS, class DH >
	std::vector< SearchInfo<MI> > SimulatedAnnealing( MS& ms, DH& dh, const int nsearch, double alpha, const int Tfactor,
																	  std::ostream& sout = std::cout, short outSI = -1, bool saveV = false,
																	  const int istart = 0 ) {
		if( alpha < 0. )
			alpha = std::pow(0.01/Tfactor,1.25/nsearch);	// emperically decided
		return SimulatedAnnealing<MI>( ms, dh, nsearch, alpha, -(double)Tfactor, sout, outSI, saveV, istart );
	}

	template < class MI, class MS, class DH >
	std::vector< SearchInfo<MI> > SimulatedAnnealing( MS& ms, DH& dh, const int nsearch, const double alpha, const double Tinit,
																	  std::ostream& sout = std::cout, short outSI = -1, bool saveV = false, 
																	  const int istart = 0 ) {
	*/
	// Simulated Allealing. Three classes required for 1) modelinfo, 2) modelspace, and 3) datahandler
	// cooltype: 0=exponential cooling, 1=linear cooling
	// outSI: -1=no-output, 0=output-all, 1=output-accepted
	// saveV: a vector of search info is returned upon request
	template < class MI, class MS, class DH >
	std::vector< SearchInfo<MI> > SimulatedAnnealing( MS& ms, DH& dh, const int nsearch, const double Tinit, 
																	  const double Tfinal, const int cooltype = 0,
																	  std::ostream& sout = std::cout, short outSI = 0, bool saveV = false,
																	  const int istart = 0 ) {
		// initialize random number generator
		_poc = 0.; float pocinc = 1./(nsearch+2);
		//std::default_random_engine generator1( std::chrono::system_clock::now().time_since_epoch().count() + std::random_device{}() );
		//std::uniform_real_distribution<float> d_uniform(0., 1.);
		//std::normal_distribution<float> d_normal(0., 1.);
		//auto rnd = std::bind( d_uniform, generator1 );
		const int nthread = omp_get_max_threads();
		Rand randA[nthread];

		// force MS, and DH to be derived from the provided interfaces at compile time
		const IModelSpace<MI>& ims = ms;
		const IDataHandler<MI>& idh = dh;

		bool outacc = outSI>=0;
		bool outrej = outSI==0;

		// cooling schedule
		const double alpha = cooltype==0 ? pow(Tfinal/Tinit, 1./(nsearch-1)) : (Tfinal-Tinit)/(nsearch-1);
		auto Cool = cooltype==0 ? 
						std::function<void(double&)>([&alpha](double &T) { T *= alpha; }) : 
						[&alpha](double &T) { T = T+alpha>0. ? T+alpha : 0.; };

		// search starts
		sout.precision(6);
		if( ! (alpha==1 && Tinit==2.0) )	// not MonteCarlo
			sout<<"### Starting simulated annealing search (#search="<<nsearch<<" alpha="<<alpha<<" Tinit="<<Tinit<<" Tfinal="<<Tfinal<<") "<<std::endl;
			sout<<"### with model state =\n"<<ms<<std::endl;

		// initial energy (force MS to be assignable to MI at compile time)
		int Ndata0;
		float E; dh.Energy( ms, E, Ndata0 );
		float Ebest = E;
		// initial temperature
		//double T = Tinit>0. ? Tinit : std::fabs(E*Tinit);
		double T = Tinit;
		// SearchInfo vector
		std::vector< SearchInfo<MI> > VSinfo; 
		VSinfo.reserve( nsearch/2*(outacc+outrej) );
		// save initial
		SearchInfo<MI> sibest( istart, T, ms, Ndata0, E, true );
		if( saveV ) VSinfo.push_back( sibest );
		sout<<sibest<<"\n";
		_poc += pocinc;

		// main loop
		const int nspawn = nthread; 
		float paccA[nspawn]; SearchInfo<MI> SIA[nspawn];
//float EaccA[nspawn];
		for( int i=istart; i<istart+nsearch; i+=nspawn ) {
			float Emin = E; int imin = -1;
		  #pragma omp parallel
		  { // spawn (simultaneously perturb) num_threads new locations from the current
			// perturb and compute Enew
			MI minew; float Enew; 
			int Ndata, isaccepted = 0;	// 0 for rejected
			try {	
				ms.Perturb( minew, false );
				dh.Energy( minew, Enew, Ndata );
				Enew *= ((float)Ndata0 / Ndata);
			} catch (const std::exception& e) {
				std::cerr<<e.what()<<"   skipped!\n";
				isaccepted = -1;	// -1 for skipped
			}
			// update Emin
			int ispawn = omp_get_thread_num();
			if( Enew < Emin ) { 
				#pragma omp critical
				if( Enew < Emin ) { Emin = Enew; imin = ispawn; }
			}
			// compute acceptance probability based on Emin
			#pragma omp barrier
//EaccA[ispawn] = Enew;
			paccA[ispawn] = isaccepted==-1 ? 0. : Paccept(Emin, Enew, T);
			// records whether the current model should be 'accepted' according to Paccept
			// note, however, that this does not decide which of the spawned models will be accepted as the new location
			if( randA[ispawn].Uniform()<paccA[ispawn] ) isaccepted = 1;
			// save searching info of current location
			SIA[ispawn] = { i+1+ispawn, T, minew, Ndata, Enew, isaccepted };
		  } // spawn ends
			// update the best
			if( Emin<Ebest ) { Ebest = Emin; sibest = SIA[imin]; }
			// include the old model as one of the candidates
			float psum = Paccept(Emin, E, T);
			// select the new location base on paccA
			// normalize paccA
			for(int ip=0; ip<nspawn; ip++) psum += paccA[ip];
			for(int ip=0; ip<nspawn; ip++) paccA[ip] /= psum;
			// select model with a random probability and accumulated pacc
			float p = randA[0].Uniform(), Paccu = 0.;
			int ip; for(ip=0; ip<nspawn; ip++) {
				Paccu += paccA[ip]; if(Paccu>p) break;
			}
			// update Energy and model state if one of the new spawns is accepted
//float Eold = E;
			if( ip != nspawn ) { const auto &si = SIA[ip]; ms.SetMState( si.info ); E = si.E; }
			// output
			if( outacc ) {
				for(const auto &si : SIA)
					if(outrej || si.accepted==1) {
						if(saveV) VSinfo.push_back( si );
						sout<<si<<"\n";
					}
//sout<<"debug info: Eold="<<Eold<<" Emin="<<Emin<<" Enew="<<E<<"(ispawn="<<ip<<") psum="<<psum<<" prand="<<p<<" (E,paccA)s: ";
//for(int i=0; i<nspawn; i++) sout<<EaccA[i]<<","<<paccA[i]<<"  "; sout<<"\n";
				if( i % 100 == 0 ) sout.flush();
			}
			// temperature decrease
			for(int it=0; it<nspawn; it++) Cool(T);
			_poc += pocinc*nspawn;	// update perc-of-completion
		}
		/*
		RWLock rwlock;
		#pragma omp parallel for schedule(dynamic, 1) shared private(Ndata)
		for( int i=istart; i<istart+nsearch; i++ ) {
			float Enew; int isaccepted = 0;	// 0 for rejected
			MI minew;
			try {	
				rwlock.lock4Read(); ms.Perturb( minew, false ); rwlock.unlock4Read();
				dh.Energy( minew, Enew, Ndata );
				Enew *= ((float)Ndata0 / Ndata);
			} catch (const std::exception& e) {
				std::cerr<<e.what()<<"   skipped!\n";
				isaccepted = -1;	// -1 for skipped
			}
			SearchInfo<MI> si( i+1, T, minew, Ndata, Enew, isaccepted );
			#pragma omp critical (MS)
			{ // critical begins
			//if( Enew<E || rnd()<exp((E-Enew)/T) ) {
			if( isaccepted==0 && rnd()<Pccept(E, Enew, T) ) {
				// update Energy and model state
				rwlock.lock4Write(); ms.SetMState( minew ); rwlock.unlock4Write();
				E = Enew;
				// mark as accepted
				si.accepted = 1;
				// update (if is) the best
				if( Enew < Ebest ) {
					Ebest = Enew;
					sibest = si;
				}
				if( outacc ) {
					// output search info
					if(saveV) VSinfo.push_back( si );
					sout<<si<<"\n";
				}
			} else if( outrej ) {
				if(saveV) VSinfo.push_back( si );
				sout<<si<<"\n";
			}
			// temperature decrease
			Cool(T);
			_poc += pocinc;	// update perc-of-completion
			if( i % 100 == 0 ) sout.flush();
			} // critical ends
		}
		*/


		// output final result
		sibest.accepted = 2;	// 2 for best fitting model
		sout<<sibest<<std::endl;
		std::sort( VSinfo.begin(), VSinfo.end() );
		if(saveV) VSinfo.push_back( sibest );
		// set model state to the best fitting model
		ms.SetMState( sibest.info );
		_poc = -1.;
		return VSinfo;
	}


	// Monte Carlo search: simulated annealing with 1) temperature fixed at 2 and 2) all search info returned
	template < class MI, class MS, class DH >
	std::vector< SearchInfo<MI> > MonteCarlo( MS& ms, DH& dh, const int nsearch, std::ostream& sout = std::cout ) {
		std::cout<<"### Monte Carlo search started. (#search = "<<nsearch<<") ###"<<std::endl;
		return SimulatedAnnealing<MI>( ms, dh, nsearch, 2.0, 2.0, 1, sout, 0 );		// save all search info
	}

	template < class MI, class MS, class DH >
	void MonteCarlo( MS& ms, DH& dh, const int nsearch, const std::string& outname ) {
		std::cout<<"### Monte Carlo search in process... (#search="<<nsearch<<") ###\n"
					<<"### results being streamed into file "<<outname<<" ###"<<std::endl;
		std::ofstream fout( outname, std::ofstream::app );
		omp_set_nested(true);
		#pragma omp parallel sections
		{	// parallel S
			#pragma omp section
			SimulatedAnnealing<MI>( ms, dh, nsearch, 2.0, 2.0, 1, fout, 0, false );	// do not save search info
			#pragma omp section
			{	// section S
			std::this_thread::sleep_for( std::chrono::seconds(1) );
			while( _poc >= 0. ) {
				std::cout<<"*** In process... "<<std::setprecision(1)<<std::setw(4)<<std::fixed<<_poc*100<<"\% completed... ***"<<std::endl<<"\x1b[A";
				std::this_thread::sleep_for( std::chrono::seconds(20) );
			}
			std::cout<<"### 100.0\% completed ###\t\t\t\n";
			}	// section E
		}	// parallel E
	}
}

#endif
