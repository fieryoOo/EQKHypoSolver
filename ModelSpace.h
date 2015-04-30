#ifndef MODELSPACE_H
#define MODELSPACE_H

#include "MyOMP.h"
#include "Searcher.h"
#include "EQKAnalyzer.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <random>
#include <chrono>
#include <thread>

class ModelSpace : public ModelInfo, public Searcher::IModelSpace<ModelInfo> {
	public:
		// construct through ModelInfo
		ModelSpace( const ModelInfo mi = ModelInfo() ) 
			: ModelInfo(mi) {
			Initialize();
		}

		// construct through param file
		ModelSpace( const std::string& fname ) {
			Initialize();
			LoadParams( fname );
		}

		// IO
		void LoadParams( const std::string& fname ) {
			std::ifstream fin(fname);
			if( ! fin )
				throw std::runtime_error("bad file");
			int nparam = 0;
			for( std::string line; std::getline(fin, line); ) {
				std::istringstream ss(line);
				if( ! (ss>>line) ) continue;
				bool succeed = false;
				if( line == "lon" ) { 
					succeed = ss >> lon; 
					if( succeed && lon<0.) lon += 360.; 
				}
				else if( line == "lat" ) succeed = ss >> lat;
				else if( line == "t0" ) succeed = ss >> t0;
				else if( line == "stk" ) succeed = ss >> stk;
				else if( line == "dip" ) succeed = ss >> dip;
				else if( line == "rak" ) succeed = ss >> rak;
				else if( line == "dep" ) succeed = ss >> dep;
				else if( line == "Rlon") succeed = ss >> Rlon;
				else if( line == "Rlat") succeed = ss >> Rlat;
				else if( line == "Rstk") succeed = ss >> Rstk;
				else if( line == "Rdip") succeed = ss >> Rdip;
				else if( line == "Rrak") succeed = ss >> Rrak;
				else if( line == "Rdep") succeed = ss >> Rdep;
				else if( line == "pertfactor") succeed = ss >> pertfactor;
				else continue;
				nparam ++;
			}
			fin.close();

			std::cout<<"### ModelSpace::LoadParams: "<<nparam<<" succed loads from param file "<<fname<<". ###"<<std::endl;

			// reset model center, perturbation ranges, and random number generators
			Initialize();
		}

		inline void SetMState( const ModelInfo& mi ) {
			static_cast<ModelInfo&>(*this) = mi;
		}

		void SetSpace( const float Clonin, const float Clatin, const float Ctimin,
				const float Cstkin, const float Cdipin, const float Crakin, const float Cdepin,
				const float Rlonin, const float Rlatin, const float Rtimin,
				const float Rstkin, const float Rdipin, const float Rrakin, const float Rdepin ) {
			Clon = Clonin; Clat = Clatin; Ctim = Ctimin;
			Cstk = Cstkin; Cdip = Cdipin; Crak = Crakin; Cdep = Cdepin;
			if( Clon<0. ) Clon+=360.;
			if( Clon<0||Clon>=360. || Clat<-90.||Clat>90. ||
					Cstk<0||Cstk>=360. || Cdip<0.||Cdip>90. || Crak<-180.||Crak>=180. || Cdep<0. )
				throw std::runtime_error("SetC");
			Rlon = Rlonin; Rlat = Rlatin; Rtim = Rtimin;
			Rstk = Rstkin; Rdip = Rdipin; Rrak = Rrakin; Rdep = Rdepin;
			if( Rlon<=0. || Rlat<=0. || Rtim<=0. || Rstk<=0. || Rdip<=0. || Rrak<=0. || Rdep<=0. )
				throw std::runtime_error("SetR");
			//validS = true;
		}

		// model space re-shape operations
		void SetPerturb( const float Plonin, const float Platin, const float Ptimin,
				const float Pstkin, const float Pdipin, const float Prakin, const float Pdepin ) {
			Plon = Plonin; Plat = Platin; Ptim = Ptimin;
			Pstk = Pstkin; Pdip = Pdipin; Prak = Prakin; Pdep = Pdepin;
			if( Plon<0. || Plat<0. || Ptim<0. || Pstk<0. || Pdip<0. || Prak<0. || Pdep<0. )
				throw std::runtime_error("Error(SetPerturb): negative purtabation!");
			if( Plon + Plat + Ptim + Pstk + Pdip + Prak + Pdep == 0 )
				throw std::runtime_error("Error(SetPerturb): model non purtabble!");
			//validP = true;
		}

		void FixEpic() {
			Plon = Plat = Ptim = 0.;
			Pstk = pertfactor * Rstk; Pdip = pertfactor * Rdip;
			Prak = pertfactor * Rrak; Pdep = pertfactor * Rdep;
		}
		void FixFocal() {
			Ptim = pertfactor * Rtim;
			Plon = pertfactor * Rlon; Plat = pertfactor * Rlat;
			Pstk = Pdip = Prak = Pdep = 0.;
		}
		void unFix() {
			Ptim = pertfactor * Rtim;
			Plon = pertfactor * Rlon; Plat = pertfactor * Rlat;
			Pstk = pertfactor * Rstk; Pdip = pertfactor * Rdip;
			Prak = pertfactor * Rrak; Pdep = pertfactor * Rdep;
		}

		void SetFreeFocal( bool rand_init = false) {
			// searching centers and ranges
			Cstk = 180.; Rstk = 180.; //Pstk = 72.;
			Cdip = 45.; Rdip = 45.; //Pdip = 18.;
			Crak = 0.; Rrak = 180.; //Prak = 72.;
			Cdep = 30.; Rdep = 30.; //Pdep = 12.;
			resetPerturb();
			// starting position
			if( rand_init ) {
				auto& rand_t = randO[omp_get_thread_num()];
				stk = rand_t.Uniform() * 360.;
				dip = rand_t.Uniform() * 90.;
				rak = rand_t.Uniform() * 360. - 180.;
				dep = rand_t.Uniform() * 60.;
			}
		}

		void BoundFocal( float Rfactor = 1. ) {
			ModelSpace ms;	// create default model space
			Rstk = ms.Rstk*Rfactor; Rdip = ms.Rdip*Rfactor;
			Rrak = ms.Rrak*Rfactor; Rdep = ms.Rdep*Rfactor;
			// reset model center, perturbation ranges, and random number generators
			Initialize();
		}

		// decide perturb step length for each parameter based on the model sensitivity to them
		// perturb steps are defined to be (ub-lb) * Sfactor, where ub and lb are the boundaries decided by:
		// assuming current model state to be the best fitting model, move away
		// from this state until the probability of acceptance <= Pthreshold
		void EstimatePerturbs( const EQKAnalyzer& eka ) {
			int Ndata;
			float Emin = eka.Energy(*this, Ndata);
			#pragma omp parallel
			{ // parallel begins
				#pragma omp sections
				{ // omp sections begins
					#pragma omp section
					{
					float lb_old = stk - Rstk, ub_old = stk + Rstk;
					float lb_stk = SearchBound( eka, stk, lb_old, Pthreshold, Emin, 10 );
					float ub_stk = SearchBound( eka, stk, ub_old, Pthreshold, Emin, 10 );
					Pstk = (ub_stk-lb_stk) * Sfactor;
					} // section 1
					#pragma omp section
					{
					float lb_old = dip - Rdip, ub_old = dip + Rdip;
					float lb_dip = SearchBound( eka, dip, lb_old, Pthreshold, Emin, 10 );
					float ub_dip = SearchBound( eka, dip, ub_old, Pthreshold, Emin, 10 );
					Pdip = (ub_dip-lb_dip) * Sfactor;
					} // section 2
					#pragma omp section
					{
					float lb_old = rak - Rrak, ub_old = rak + Rrak;
					float lb_rak = SearchBound( eka, rak, lb_old, Pthreshold, Emin, 10 );
					float ub_rak = SearchBound( eka, rak, ub_old, Pthreshold, Emin, 10 );
					Prak = (ub_rak-lb_rak) * Sfactor;
					} // section 3
					#pragma omp section
					{
					float lb_old = dep - Rdep, ub_old = dep + Rdep;
					float lb_dep = SearchBound( eka, dep, lb_old, Pthreshold, Emin, 10 );
					float ub_dep = SearchBound( eka, dep, ub_old, Pthreshold, Emin, 10 );
					Pdep = (ub_dep-lb_dep) * Sfactor;
					} // section 4
					#pragma omp section
					{
					float lb_old = lon - Rlon, ub_old = lon + Rlon;
					float lb_lon = SearchBound( eka, lon, lb_old, Pthreshold, Emin, 10 );
					float ub_lon = SearchBound( eka, lon, ub_old, Pthreshold, Emin, 10 );
					Plon = (ub_lon-lb_lon) * Sfactor;
					} // section 5
					#pragma omp section
					{
					float lb_old = lat - Rlat, ub_old = lat + Rlat;
					float lb_lat = SearchBound( eka, lat, lb_old, Pthreshold, Emin, 10 );
					float ub_lat = SearchBound( eka, lat, ub_old, Pthreshold, Emin, 10 );
					Plat = (ub_lat-lb_lat) * Sfactor;
					} // section 6
					#pragma omp section
					{
					float lb_old = t0 - Rtim, ub_old = t0 + Rtim;
					float lb_tim = SearchBound( eka, t0, lb_old, Pthreshold, Emin, 10 );
					float ub_tim = SearchBound( eka, t0, ub_old, Pthreshold, Emin, 10 );
					Ptim = (ub_tim-lb_tim) * Sfactor;
					} // section 7
				} // omp sections ends
			} // parallel ends
			std::cout<<"*** Model state after estemating perturb step sizes: "<<*this<<std::endl;
		}

		void Perturb( ModelInfo& minew ) const {
			//if( ! (validS && validP) ) throw std::runtime_error("incomplete model space");
			if( ! isValid() ) throw std::runtime_error("incomplete model info");

			// stk
			float stk_cur = ShiftInto(this->stk, Cstk-Rstk, Cstk+Rstk, 360.);
			if( Rstk >= 180. ) {
				minew.stk = Neighbour_Cycle(stk_cur, Pstk, 0., 360.);
			} else {
				minew.stk = Neighbour_Reflect(stk_cur, Pstk, Cstk-Rstk, Cstk+Rstk);
				minew.stk = ShiftInto(minew.stk, 0., 360., 360.);
			}

			// dip
			float lb = Cdip-Rdip, ub = Cdip+Rdip;
			if( lb < 0. ) lb = 0.;
			if( ub > 90. ) ub = 90.;
			minew.dip = Neighbour_Reflect(this->dip, Pdip, lb, ub);

			// rak
			float rak_cur = ShiftInto(this->rak, Crak-Rrak, Crak+Rrak, 360.);
			if( Rrak >= 180. ) {
				minew.rak = Neighbour_Cycle(rak_cur, Prak, -180., 180.);
			} else {
				minew.rak = Neighbour_Reflect(rak_cur, Prak, Crak-Rrak, Crak+Rrak);
				minew.rak = ShiftInto(minew.rak, -180., 180., 360.);
			}

			// dep
			lb = Cdep-Rdep; ub = Cdep+Rdep;
			if( lb < 0. ) lb = 0.;
			minew.dep = Neighbour_Reflect(this->dep, Pdep, lb, ub);

			// longitude
			minew.lon = Neighbour_Reflect(this->lon, Plon, Clon-Rlon, Clon+Rlon);

			// latitude
			minew.lat = Neighbour_Reflect(this->lat, Plat, Clat-Rlat, Clat+Rlat);

			// origin time
			minew.t0 = Neighbour_Reflect(this->t0, Ptim, Ctim-Rtim, Ctim+Rtim);
		}

		/* streaming perturbation ranges */
		friend std::ostream& operator<< ( std::ostream& o, ModelSpace& ms ) {
			o<<"  lon ("<<ms.Clon-ms.Rlon<<"~"<<ms.Clon+ms.Rlon<<", "<<ms.Plon<<")"
			<<"  lat ("<<ms.Clat-ms.Rlat<<"~"<<ms.Clat+ms.Rlat<<", "<<ms.Plat<<")"
			<<"  tim ("<<ms.Ctim-ms.Rtim<<"~"<<ms.Ctim+ms.Rtim<<", "<<ms.Ptim<<")\n"
			<<"  stk ("<<ms.Cstk-ms.Rstk<<"~"<<ms.Cstk+ms.Rstk<<", "<<ms.Pstk<<")"
			<<"  dip ("<<ms.Cdip-ms.Rdip<<"~"<<ms.Cdip+ms.Rdip<<", "<<ms.Pdip<<")"
			<<"  rak ("<<ms.Crak-ms.Rrak<<"~"<<ms.Crak+ms.Rrak<<", "<<ms.Prak<<")"
			<<"  dep ("<<ms.Cdep-ms.Rdep<<"~"<<ms.Cdep+ms.Rdep<<", "<<ms.Pdep<<")";
			return o;
		}

	protected:
		static constexpr float Pthreshold = 0.005;   /* the threshold for probability in searching for parameter
																		sensitivity prior to the Monte Carlo search, */
	   static constexpr float Sfactor = 0.1;        /* and the step half-length for the search as a fraction
																		of (ub-lb) decided by Pthreshold */

	private: // variables
		//bool validS{false}, validP{false};
		float Clon{NaN}, Clat{NaN}, Ctim{NaN}, Cstk{NaN}, Cdip{NaN}, Crak{NaN}, Cdep{NaN}; // model center
		float Rlon{0.3}, Rlat{0.3}, Rtim{5.}, Rstk{30.}, Rdip{15.}, Rrak{30.}, Rdep{7.5}; // model param radius
		float Plon{NaN}, Plat{NaN}, Ptim{NaN}, Pstk{NaN}, Pdip{NaN}, Prak{NaN}, Pdep{NaN}; // perturb length ( gaussian half length )
		float pertfactor{0.1};
		mutable std::vector<Rand> randO;

	private:	// methods
		// initialize perturbation length and random number generators
		void resetPerturb() {
			// re-compute perturbation ranges
			Ptim = pertfactor * Rtim;
			Plon = pertfactor * Rlon; Plat = pertfactor * Rlat;
			Pstk = pertfactor * Rstk; Pdip = pertfactor * Rdip;
			Prak = pertfactor * Rrak; Pdep = pertfactor * Rdep;
		}
		void Initialize() {
			randO.clear();
			// produce Rand object for each thread
			for(int i=0; i<omp_get_max_threads(); i++) {
				// apply separated seed by sleeping
				randO.push_back( Rand() );
				std::this_thread::sleep_for(std::chrono::milliseconds(100));
			}
			// set model center to the current MState
			Clon = lon; Clat = lat; Ctim = t0;
			Cstk = stk; Cdip = dip; Crak = rak; Cdep = dep;
			// compute perturbation ranges
			resetPerturb();
		}

		// perturbing values
		inline float Neighbour_Cycle( float valold, float hlen, float lb, float ub ) const {
			float valnew = valold + randO[omp_get_thread_num()].Normal()*hlen;
			float range = ub - lb;
			while( valnew >= ub ) { valnew -= range; }
			while( valnew < lb ) { valnew += range; }
			return valnew;
		}
		inline float Neighbour_Reflect( float valold, float hlen, float lb, float ub ) const {
			if( hlen == 0. ) return valold;
			if( valold<lb || valold>ub ) throw std::runtime_error("old val out of boundary");
			float range = ub - lb;
			float shift = randO[omp_get_thread_num()].Normal() * hlen;
			if( shift > range ) shift = range;
			else if( shift < -range ) shift = -range;
			float valnew = valold + shift;
			if( valnew >= ub ) { valnew = 2.*ub-valnew; }
			else if( valnew < lb ) { valnew = 2.*lb-valnew; }
			return valnew;
		}

		// shift by T multiples according to lower and upper bound. Results not guranteed to be in the range
		inline float ShiftInto(float val, float lb, float ub, float T) const {
			while(val >= ub) val -= T;
			while(val < lb) val += T;
			return val;
		}

		float SearchBound( const EQKAnalyzer& eka, float& key, float bound, float Pthsd, float Emin, int nsearch ) {
			const float key_orig = key; // save the original key
			float direct = 1, steplen = (bound-key) * 0.5;
			for(int isearch=0; isearch<nsearch; isearch++) {
				key += steplen * direct;
				int Ndata;
				float E = eka.Energy(*this, Ndata);
				float P = exp(0.5*(Emin-E));
				if( P < Pthsd ) direct = -1;    // moving backward
				else direct = 1;
				steplen *= 0.5;
				//std::cerr<<isearch<<"   "<<E<<" "<<P<<"   "<<direct<<" "<<steplen<<std::endl;
			}
			bound = key + steplen * direct;
			key = key_orig;   // set the key back to its original value
			return bound;    // return final position after nsearches
		}
};


#endif