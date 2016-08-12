#ifndef MODELINFO_H
#define MODELINFO_H

#include "MyOMP.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <random>
#include <chrono>
#include <thread>

/* -------------------- the RNG class-------------------- */
/*
class Rand {
   std::default_random_engine generator1;
   std::uniform_real_distribution<float> d_uniform;
   std::normal_distribution<float> d_normal;
public:
   Rand() // add a true random number from std::random_device to the time seed to ensure thread safety 
      : generator1( std::chrono::system_clock::now().time_since_epoch().count() + std::random_device{}() )
      , d_uniform(0., 1.)
      , d_normal(0., 1.) {}

   //~Rand() {} // should not be defined!!!

   inline float Uniform() { return d_uniform(generator1); }
   inline float Normal() { return d_normal(generator1); }

};
*/

template< class T >
struct FocalInfo {
	static constexpr float NaN = -12345.;
	static constexpr float DEPMAX = 60.;
   T stk, dip, rak, dep; mutable T M0;

   //FocalInfo( T stkin = 180, T dipin = 45, T rakin = 0, T depin = 10 )
   FocalInfo( T stkin = NaN, T dipin = NaN, T rakin = NaN, T depin = NaN, T M0in = NaN )
      : stk(stkin), dip(dipin), rak(rakin), dep(depin), M0(M0in) {}

	/* check validation */
	virtual bool isValid() const {
		return (stk>=0.&&stk<360.) && (dip>=0.&&dip<=90.) &&
				 (rak>=-180.&&rak<180.) && (dep>=0.&&dep<DEPMAX) ;
	}

	/* correct into the right range */
	inline void CorrectRange() {
		Correct2PI();
		if( stk==NaN || rak==NaN || dip==NaN || dip>=180. || dip<=-90. )
			throw std::runtime_error("Error(FocalInfo::CorrectRange): Invalid finfo!");
		if( dip >= 90. ) dip = 180. - dip;
		if( dip < 0. ) dip = -dip;
	}

	/* Transfer to the auxiliary nodal plane and slip of the current */
	void Auxiliary() {
		float s2, d2, r2;
		// convert deg to rad
		float deg2rad = M_PI/180.;
		stk *= deg2rad; dip *= deg2rad; rak *= deg2rad;
		// dip
		d2 = acos( sin(rak)*sin(dip) );
		// rak
		float sin_r2 = cos(dip) / sin(d2);
		float sin_dphi = cos(rak) / sin(d2);
		float cos_r2 = -sin(dip) * sin_dphi;	// 0 ~ pi
		r2 = asin(sin_r2);	// -pi/2 ~ pi/2
		if( cos_r2 < 0. ) r2 = M_PI - r2;
		// stk
		float cos_dphi = -1. / (tan(dip)*tan(d2));
		float dphi = asin(sin_dphi);
		if( cos_dphi < 0. ) dphi = M_PI - dphi;
		s2 = stk - dphi;
		// check and correct quadrant
		if( d2 > 0.5*M_PI ) {
			s2 += M_PI;
			d2 = M_PI - d2;
			r2 = 2.*M_PI - r2;
		}
		// convert back;
		float rad2deg = 1./deg2rad;
		stk = s2*rad2deg; dip = d2*rad2deg; rak = r2*rad2deg;
		Correct2PI();
	}

	// NED(default)=NorthEastDown, USE=UpSouthEast
	void printMomentTensor( const float M0, const bool USE=false ) const {
		auto MT = MomentTensor(M0);
		printM(std::cout, MT[0], MT[1], MT[2], MT[3], MT[4], MT[5], USE)<<std::endl;
	}
	std::array<float, 6> MomentTensor( const float M0 ) const {
		float deg2rad = M_PI/180.;
		float stk = this->stk, dip = this->dip, rak = this->rak;
		stk *= deg2rad; dip *= deg2rad; rak *= deg2rad;
		float sins = sin(stk), coss = cos(stk), sin2s = sin(2.*stk), cos2s = cos(2.*stk);
		float sind = sin(dip), cosd = cos(dip), sin2d = sin(2.*dip), cos2d = cos(2.*dip);
		float sinr = sin(rak), cosr = cos(rak);
		std::array<float, 6> MT;
		MT[0] = -M0 * (sind*cosr*sin2s + sin2d*sinr*sins*sins);	// xx
		MT[1] =  M0 * (sind*cosr*sin2s - sin2d*sinr*coss*coss);	// yy
		MT[2] =  M0 * (sin2d*sinr);										// zz
		MT[3] =  M0 * (sind*cosr*cos2s + sin2d*sinr*sins*coss);	// xy
		MT[4] = -M0 * (cosd*cosr*coss + cos2d*sinr*sins);			// xz
		MT[5] = -M0 * (cosd*cosr*sins - cos2d*sinr*coss);			// yz
		return MT;
	}

	void printMTDecomposed( const float M0, const bool USE=false ) const {
		float deg2rad = M_PI/180.;
		float stk = this->stk, dip = this->dip, rak = this->rak;
		stk *= deg2rad; dip *= deg2rad; rak *= deg2rad;
		float sins = sin(stk), coss = cos(stk), sin2s = sin(2.*stk), cos2s = cos(2.*stk);
		float sind = sin(dip), cosd = cos(dip), sin2d = sin(2.*dip), cos2d = cos(2.*dip);
		float sinr = sin(rak), cosr = cos(rak);
		std::cout<<cosd*cosr<<" x\n"; printM(std::cout, 0, 0, 0, 0, -coss, -sins, USE)<<"\n";
		std::cout<<sind*cosr<<" x\n"; printM(std::cout, -sin2s, sin2s, 0, cos2s, 0, 0, USE)<<"\n";
		std::cout<<-cos2d*sinr<<" x\n"; printM(std::cout, 0, 0, 0, 0, sins, -coss, USE)<<"\n";
		std::cout<<sin2d*sinr<<" x\n"; printM(std::cout, -sins*sins, -coss*coss, 1, 0.5*sin2s, 0, 0, USE)<<std::endl;
	}

	std::ostream &printM( std::ostream &o, float M11, float M22, float M33, 
								 float M12, float M13, float M23, const bool USE=false) const {
		if( USE ) {
			float Mtmp = M11; M11 = M33; M33 = M22; M22 = Mtmp;
			Mtmp=M12; M12 = M13; M13 = -M23; M23 = -Mtmp;
		}
		o<<std::setprecision(5)<<std::fixed
		 <<std::setw(9)<<M11<<std::setw(9)<<M12<<std::setw(9)<<M13<<"\n"
		 <<std::setw(9)<<" "<<std::setw(9)<<M22<<std::setw(9)<<M23<<"\n"
		 <<std::setw(9)<<" "<<std::setw(9)<<" "<<std::setw(9)<<M33;//<<std::endl;
		return o;
	}

   friend bool operator== ( const FocalInfo<T>& fi1, const FocalInfo<T>& fi2 ) {
      T dis_st = fabs(fi1.stk - fi2.stk);
      T dis_di = fabs(fi1.dip - fi2.dip);
      T dis_ra = fabs(fi1.rak - fi2.rak);
      T dis_de = fabs(fi1.dep - fi2.dep);
      return (dis_st<toler && dis_di<toler && dis_ra<toler && dis_de<0.1*toler);
   }

   friend std::ostream &operator<< ( std::ostream& o, const FocalInfo& f ) {
      //o<<std::fixed<<std::setprecision(2)<<f.stk<<" "<<std::setw(6)<<f.dip<<" "<<std::setw(6)<<f.rak<<"  "<<std::setw(6)<<f.dep;
      o<<std::fixed<<std::setprecision(3)
		 <<std::setw(7)<<f.stk<<" "<<std::setw(7)<<f.dip<<" "<<std::setw(8)<<f.rak<<" "<<std::setw(6)<<f.dep<<" "<<std::setw(7)<<std::scientific<<f.M0; 
      return o; 
   }

   friend std::istream &operator>> ( std::istream& i, FocalInfo& f ) {
		i >> f.stk >> f.dip >> f.rak >> f.dep >> f.M0; return i;
   }

protected:
	static constexpr float toler = 1.0e-3;
private:
	/* shift by 2PIs into the correct range */
	void Correct2PI() {
		if( stk==NaN || rak==NaN ) return;
		while( stk >= 360. ) stk -= 360.;
		while( stk < 0. ) stk += 360.;
		while( rak >= 180. ) rak -= 360.;
		while( rak < -180. ) rak += 360.;
	}

};
typedef float ftype;

/* earthquake epicenter information */
struct EpicInfo {
	static constexpr float NaN = FocalInfo<ftype>::NaN;
	float lon, lat, t0;

	EpicInfo( float lonin = NaN, float latin = NaN, float t0in = NaN )
		: lon(lonin), lat(latin), t0(t0in) {}

	virtual bool isValid() const {
		return ( (lon<360.&&lon>=0.) && (lat>=-90.&&lat<=90.) && t0!=NaN );
	}

	friend bool operator== ( const EpicInfo& ei1, const EpicInfo& ei2 ) {
		float dis_lon = fabs(ei1.lon - ei2.lon);
		float dis_lat = fabs(ei1.lat - ei2.lat);
		float dis_t = fabs(ei1.t0 - ei2.t0);
		return (dis_lon<toler && dis_lat<toler && dis_t<10.*toler);
	}

	friend std::ostream& operator<< ( std::ostream& o, const EpicInfo &e ) {
		o<<std::fixed<<std::setprecision(4)<<e.lon<<" "<<e.lat<<" "<<std::setw(7)<<e.t0; 
		return o; 
	}

   friend std::istream &operator>> ( std::istream& i, EpicInfo &e ) {
		i >> e.lon >> e.lat >> e.t0; return i;
   }

protected:
	static constexpr float toler = 1.0e-5;
};


/* model space */
class ModelInfo : public FocalInfo<ftype>, public EpicInfo {
	public:
		using FocalInfo::NaN;

		ModelInfo() {}

		ModelInfo( const float lonin, const float latin, const float timin,
				const float stkin, const float dipin, const float rakin, const float depin, const float M0in = 1. )
			: EpicInfo(lonin, latin, timin)
			, FocalInfo(stkin, dipin, rakin, depin, M0in) {}

		ModelInfo( const std::string& line ) {
			Loadline(line);
		}

		ModelInfo( const EpicInfo& einfo, const FocalInfo& finfo )
			: EpicInfo(einfo), FocalInfo(finfo) {}

		void Loadline( const std::string& line ) {
			std::stringstream ss(line);
			if( !( ss >> lon >> lat >> t0 >> stk >> dip >> rak >> dep ) )
				throw std::runtime_error("Error(ModelInfo::ModelInfo): format error in: " + line);
			ss >> M0;
		}

		virtual bool isValid() const {
			return ( FocalInfo<ftype>::isValid() && EpicInfo::isValid() );
		}

		virtual void Correct() {
			float stkold=stk, rakold=rak;
			stk = ShiftInto( stk, 0., 360., 360. );	
			rak = ShiftInto( rak, -180., 180., 360. );
			if( ! isValid() )
				throw std::runtime_error("Error(ModelInfo::Correct): invalid model info = "+toString());
			if( stk!=stkold || rak!=rakold )
				std::cerr<<"Warning(ModelInfo::Correct): either stk or rak got corrected ("<<*this<<")!"<<std::endl;
			//dip = BoundInto( dip, 0., 90. );
			//dep = BoundInto( dep, 0., DEPMAX );
		}

		std::string toString() const {
			std::stringstream ss; ss << *this;
			return ss.str();
		}

		friend bool operator== ( const ModelInfo& ms1, const ModelInfo& ms2 ) {
			return ( ( static_cast< const FocalInfo<ftype>& >(ms1) == static_cast< const FocalInfo<ftype>& >(ms2) )
					&& ( static_cast< const EpicInfo& >(ms1) == static_cast< const EpicInfo& >(ms2) ) );
		}

		friend std::ostream &operator<< ( std::ostream& o, const ModelInfo& m ) {
			o<< static_cast<const FocalInfo<ftype>&>(m) << "   " << static_cast<const EpicInfo&>(m);
			return o;
		}

		friend std::istream &operator>> ( std::istream& i, ModelInfo& m ) {
			i >> static_cast<FocalInfo<ftype>&>(m) >> static_cast<EpicInfo&>(m); return i;
		}

	private:
		// shift by T multiples according to lower and upper bound. Results not guranteed to be in the range
		inline float ShiftInto( float val, float lb, float ub, float T) const {
			while(val >= ub) val -= T;
			while(val < lb) val += T;
			return val;
		}
		inline float BoundInto( float val, float lb, float ub ) const {
			return val<lb ? lb : (val>ub ? ub : val);
		}
};

#endif
