#ifndef MODELSPACE_H
#define MODELSPACE_H

//#include "FocalSearcher.h"
#include "MyOMP.h"
#include <chrono>
#include <thread>


/* -------------------- the RNG class-------------------- */
class Rand {
   std::default_random_engine generator1;
   std::uniform_real_distribution<float> d_uniform;
   std::normal_distribution<float> d_normal;
public:
   Rand() /* add a true random number from std::random_device to the time seed to ensure thread safety */
      : generator1( std::chrono::system_clock::now().time_since_epoch().count() + std::random_device{}() )
      , d_uniform(0., 1.)
      , d_normal(0., 1.) {}

   //~Rand() {} // should not be defined!!!

   inline float Uniform() { return d_uniform(generator1); }
   inline float Normal() { return d_normal(generator1); }

};

template< class T >
struct FocalInfo {
	static constexpr float NaN = -12345.;
   T strike, dip, rake, depth;

   //FocalInfo( T strikein = 180, T dipin = 45, T rakein = 0, T depthin = 10 )
   FocalInfo( T strikein = NaN, T dipin = NaN, T rakein = NaN, T depthin = NaN )
      : strike(strikein), dip(dipin), rake(rakein), depth(depthin) {}

	/* check validation */
	virtual bool isValid() const {
		return (strike>=0.&&strike<360.) && (dip>=0.&&dip<=90.) &&
				 (rake>=-180.&&rake<180.) && (depth>=0.) ;
	}

	/* correct into the right range */
	inline void CorrectRange() {
		Correct2PI();
		if( strike==NaN || rake==NaN || dip==NaN || dip>=180. || dip<=-90. )
			throw std::runtime_error("Error(FocalInfo::CorrectRange): Invalid finfo!");
		if( dip >= 90. ) dip = 180. - dip;
		if( dip < 0. ) dip = -dip;
	}

	/* Transfer to the auxiliary nodal plane and slip of the current */
	void Auxiliary() {
		float s2, d2, r2;
		// convert deg to rad
		float deg2rad = M_PI/180.;
		strike *= deg2rad; dip *= deg2rad; rake *= deg2rad;
		// dip
		d2 = acos( sin(rake)*sin(dip) );
		// rake
		float sin_r2 = cos(dip) / sin(d2);
		float sin_dphi = cos(rake) / sin(d2);
		float cos_r2 = -sin(dip) * sin_dphi;	// 0 ~ pi
		r2 = asin(sin_r2);	// -pi/2 ~ pi/2
		if( cos_r2 < 0. ) r2 = M_PI - r2;
		// strike
		float cos_dphi = -1. / (tan(dip)*tan(d2));
		float dphi = asin(sin_dphi);
		if( cos_dphi < 0. ) dphi = M_PI - dphi;
		s2 = strike - dphi;
		// check and correct quadrant
		if( d2 > 0.5*M_PI ) {
			s2 += M_PI;
			d2 = M_PI - d2;
			r2 = 2.*M_PI - r2;
		}
		// convert back;
		float rad2deg = 1./deg2rad;
		strike = s2*rad2deg; dip = d2*rad2deg; rake = r2*rad2deg;
		Correct2PI();
	}

   friend std::ostream& operator<< ( std::ostream& o, const FocalInfo& f ) {
      o.precision(6);
      o<<std::setw(6)<<f.strike<<" "<<std::setw(6)<<f.dip<<" "<<std::setw(6)<<f.rake<<"  "<<std::setw(6)<<f.depth; 
      return o; 
   }

   friend bool operator== ( const FocalInfo<T>& fi1, const FocalInfo<T>& fi2 ) {
      T dis_st = fabs(fi1.strike - fi2.strike);
      T dis_di = fabs(fi1.dip - fi2.dip);
      T dis_ra = fabs(fi1.rake - fi2.rake);
      T dis_de = fabs(fi1.depth - fi2.depth);
      return (dis_st<0.1 && dis_di<0.1 && dis_ra<0.1 && dis_de<0.1);
   }

private:
	/* shift by 2PIs into the correct range */
	void Correct2PI() {
		if( strike==NaN || rake==NaN ) return;
		while( strike >= 360. ) strike -= 360.;
		while( strike < 0. ) strike += 360.;
		while( rake >= 180. ) rake -= 360.;
		while( rake < -180. ) rake += 360.;
	}

};
typedef float ftype;

/* earthquake epicenter information */
struct EpicInfo {
	static constexpr float NaN = FocalInfo<ftype>::NaN;
	float lon, lat, t0;

	EpicInfo( float lonin = -12345., float latin = -12345., float t0in = -12345. )
		: lon(lonin), lat(latin), t0(t0in) {}

	virtual bool isValid() const {
		return ( (lon<360.&&lon>=0.) && (lat>=-90.&&lat<=90.) && t0!=NaN );
	}

	friend std::ostream& operator<< ( std::ostream& o, const EpicInfo& e ) {
		o.precision(6);
		o<<std::setw(6)<<e.lon<<" "<<std::setw(6)<<e.lat<<"  "<<std::setw(6)<<e.t0; 
		return o; 
	}

	friend bool operator== ( const EpicInfo& ei1, const EpicInfo& ei2 ) {
		float dis_lon = fabs(ei1.lon - ei2.lon) * 100.;
		float dis_lat = fabs(ei1.lat - ei2.lat) * 100.;
		float dis_t = fabs(ei1.t0 - ei2.t0);
		return (dis_lon<0.01 && dis_lat<0.01 && dis_t<0.01);
	}
};


/* model space */
class ModelInfo : public FocalInfo<ftype>, public EpicInfo {
	public:
		using FocalInfo::NaN;

		ModelInfo() {}

		ModelInfo( const float lonin, const float latin, const float timin,
				const float stkin, const float dipin, const float rakin, const float depin )
			: EpicInfo(lonin, latin, timin)
			  , FocalInfo(stkin, dipin, rakin, depin) {}

		ModelInfo( const EpicInfo& einfo, const FocalInfo& finfo )
			: EpicInfo(einfo), FocalInfo(finfo) {}

		virtual bool isValid() const {
			return ( FocalInfo<ftype>::isValid() && EpicInfo::isValid() );
		}

		friend std::ostream& operator<< ( std::ostream& o, const ModelInfo& m ) {
			o<< static_cast< const FocalInfo<ftype>& >(m) << "   " << static_cast< const EpicInfo& >(m);
			return o;
		}

		friend bool operator== ( const ModelInfo& ms1, const ModelInfo& ms2 ) {
			return ( ( static_cast< const FocalInfo<ftype>& >(ms1) == static_cast< const FocalInfo<ftype>& >(ms2) )
					&& ( static_cast< const EpicInfo& >(ms1) == static_cast< const EpicInfo& >(ms2) ) );
		}
};

class ModelSpace : public ModelInfo {
	public:
		ModelSpace()
			: validS(false), validP(false)
			  , Rlon(NaN), Rlat(NaN), Rtim(NaN)
			  , Plon(NaN), Plat(NaN), Ptim(NaN)
			  , Rstk(NaN), Rdip(NaN), Rrak(NaN), Rdep(NaN)
			  , Pstk(NaN), Pdip(NaN), Prak(NaN), Pdep(NaN) {
				  // produce Rand object for each thread
				  for(int i=0; i<omp_get_max_threads(); i++) {
					  // apply separated seed by sleeping
					  randO.push_back( Rand() );
					  std::this_thread::sleep_for(std::chrono::milliseconds(100));
				  }
			  }

		ModelSpace( const ModelInfo& mi ) 
			: ModelSpace() {
				strike = mi.strike;
				dip = mi.dip;
				rake = mi.rake;
				depth = mi.depth;
				lon = mi.lon;
				lat = mi.lat;
				t0 = mi.t0;
			}

		inline void SetMState( const ModelInfo& mi ) {
			*( static_cast< ModelInfo* >(this) ) = mi;
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
			validS = true;
		}

		void SetPerturb( const float Plonin, const float Platin, const float Ptimin,
				const float Pstkin, const float Pdipin, const float Prakin, const float Pdepin ) {
			Plon = Plonin; Plat = Platin; Ptim = Ptimin;
			Pstk = Pstkin; Pdip = Pdipin; Prak = Prakin; Pdep = Pdepin;
			if( Plon<=0. || Plat<=0. || Ptim<=0. || Pstk<=0. || Pdip<=0. || Prak<=0. || Pdep<=0. )
				throw std::runtime_error("Error(SetPerturb): negative purtabation!");
			validP = true;
		}

		void SetFreeFocal() {
			// searching centers and ranges
			Cstk = 180.; Rstk = 180.; //Pstk = 72.;
			Cdip = 45.; Rdip = 45.; //Pdip = 18.;
			Crak = 0.; Rrak = 180.; //Prak = 72.;
			Cdep = 30.; Rdep = 30.; //Pdep = 12.;
			// starting position
			auto& rand_t = randO[omp_get_thread_num()];
			strike = rand_t.Uniform() * 360.;
			dip = rand_t.Uniform() * 90.;
			rake = rand_t.Uniform() * 360. - 180.;
			depth = rand_t.Uniform() * 60.;
		}

		inline float Neighbour_Reflect( float valold, float hlen, float lb, float ub ) {
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
		inline float ShiftInto(float val, float lb, float ub, float T) {
			while(val >= ub) val -= T;
			while(val < lb) val += T;
			return val;
		}

		void Perturb( ModelInfo& minew ) {
			if( ! (validS && validP) ) throw std::runtime_error("invalid");

			float strike_cur = ShiftInto(this->strike, Cstk-Rstk, Cstk+Rstk, 360.);
			minew.strike = Neighbour_Reflect(strike_cur, Pstk, Cstk-Rstk, Cstk+Rstk);
			minew.strike = ShiftInto(minew.strike, 0., 360., 360.);

			float lb = Cdip-Rdip, ub = Cdip+Rdip;
			if( lb < 0. ) lb = 0.;
			if( ub > 90. ) ub = 90.;
			minew.dip = Neighbour_Reflect(this->dip, Pdip, lb, ub);

			float rake_cur = ShiftInto(this->rake, Crak-Rrak, Crak+Rrak, 360.);
			minew.rake = Neighbour_Reflect(rake_cur, Prak, Crak-Rrak, Crak+Rrak);
			minew.rake = ShiftInto(minew.rake, -180., 180., 360.);

			lb = Cdep-Rdep; ub = Cdep+Rdep;
			if( lb < 0. ) lb = 0.;
			minew.depth = Neighbour_Reflect(this->depth, Pdep, lb, ub);

			minew.lon = Neighbour_Reflect(this->lon, Plon, Clon-Rlon, Clon+Rlon);

			minew.lat = Neighbour_Reflect(this->lat, Plat, Clat-Rlat, Clat+Rlat);

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
	}

	private:
		bool validS, validP;
		float Clon, Clat, Ctim, Cstk, Cdip, Crak, Cdep; // model center
		float Rlon, Rlat, Rtim, Rstk, Rdip, Rrak, Rdep; // model param radius
		float Plon, Plat, Ptim, Pstk, Pdip, Prak, Pdep; // perturb length ( gaussian half length )
		std::vector<Rand> randO;
};


template < class T >
struct SearchInfo {
	int isearch;
	float E;
	T info;
};

#endif
