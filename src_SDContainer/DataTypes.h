#ifndef DATATYPES_H
#define DATATYPES_H

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

/* data structure for group, phase, and amplitude data and their std-dev
	at a single azimuth */
struct AziData {
   //bool valid = false;
   float azi = NaN;
	float Gdata = NaN, Pdata = NaN, Adata = NaN;
	float user = NaN;	// one extra field for amplitude
	//float Gstd = NaN, Pstd = NaN, Astd = NaN;

   AziData( float ini = NaN ) 
		: azi(ini), Gdata(ini), Pdata(ini), Adata(ini), user(ini) {}

   AziData( float aziin, float Gdatain, float Pdatain, float ampin, float userin=NaN )
      : azi(aziin), Gdata(Gdatain), Pdata(Pdatain), Adata(ampin), user(userin) {} //, valid(validin)

   //bool isGood() const { return ( valid && Gdata!=NaN && Gstd!=NaN && Pdata!=NaN && Pstd!=NaN && Adata!=NaN ); }


	// check if (all attributes are) within a given range defined by two other AziDatas
	virtual bool isWithin( const AziData& lbound, const AziData& ubound ) const {
		return (
			azi >= lbound.azi && azi <= ubound.azi &&
			Gdata >= lbound.Gdata && Gdata <= ubound.Gdata &&
			Pdata >= lbound.Pdata && Pdata <= ubound.Pdata &&
			Adata >= lbound.Adata && Adata <= ubound.Adata
		);
	}
	//if(debug) std::cerr<<"ExcludeBad: "<<adcur.azi<<"   "<<adcur.misG<<" "<<lbG<<" "<<ubG<<"   "<<adcur.misP<<" "<<lbP<<" "<<ubP<<"   "<<adcur.A<<" "<<lbA<<" "<<ubA<<std::endl;
	/* friends */
	/*
	friend void Suppress2PI( std::vector<AziData>& data, float per ) {
		if( data.size() <= 1 ) return;
      float pero2 = per * 0.5;
		// find mediam of Pdata
      std::vector<float> datatmp;
		for( const auto& ad : data )
			datatmp.push_back(ad.Pdata);
      int halfs = datatmp.size() / 2;
      std::nth_element(datatmp.begin(), datatmp.begin()+halfs-1, datatmp.end());
      float median = datatmp.at(halfs-1);
		datatmp.clear();
		// correct 2pis
		float lbound = median - pero2, ubound = median + pero2;
		for( auto& ad : data ) {
			float& adPdata = ad.Pdata;
         while( adPdata > ubound ) adPdata -= per;
         while( adPdata <= lbound ) adPdata += per;
      }

	} */

	/* math operators */

   friend bool operator<( const AziData& ad1, const AziData& ad2 ) {
      return (ad1.azi < ad2.azi);
   }

	// unary negation 
	AziData operator-() const { return AziData( -azi, -Gdata, -Pdata, -Adata, -user ); }

	// addition (Be careful! azimuths get summed too!)
	AziData& operator+=( const AziData& ad2 ) {
		azi += ad2.azi; Gdata += ad2.Gdata;
		Pdata += ad2.Pdata; Adata += ad2.Adata;
		user += ad2.user;
		return *this; 
	}
	friend AziData operator+( const AziData& ad1, const AziData& ad2 ) {
		AziData adres = ad1;
		adres += ad2;
		return adres;
	}

	// subtraction
   AziData& operator-=( const AziData& ad2 ) {
		(*this) += -ad2;
      return *this;
   }
   friend AziData operator-( const AziData& ad1, const AziData& ad2 ) {
      AziData adres = ad1;
      adres -= ad2;
      return adres;
   }

	// multiplication (*float)
	AziData& operator*=( const float mul ) { 
		azi *= mul; Gdata *= mul; Pdata *= mul; Adata *= mul; user *= mul;
		return *this; 
	}
	friend AziData operator*( const AziData& ad1, float mul ) {
		AziData adres = ad1;
		adres *= mul;
		return adres;
	}
	friend AziData operator*( float mul, const AziData& ad1 ) {
		AziData adres = ad1;
		adres *= mul;
		return adres;
	}

	// multiplication (*AziData)
	AziData& operator*=( const AziData& ad2 ) { 
		azi *= ad2.azi; Gdata *= ad2.Gdata; 
		Pdata *= ad2.Pdata; Adata *= ad2.Adata;
		user *= ad2.user;
		return *this; 
	}
	friend AziData operator*( const AziData& ad1, const AziData& ad2 ) {
		AziData adres = ad1;
		adres *= ad2;
		return adres;
	}

	// division (*float)
	AziData& operator/=( const float den ) {
		float mul = 1./den;
		*this *= mul;
		return *this;
	}
	friend AziData operator/( const AziData& ad1, float den ) {
		AziData adres = ad1;
		adres /= den;
		return adres;
	}
	friend AziData operator/( float den, const AziData& ad1 ) {
		AziData adres;
		adres.azi = den / ad1.azi;
		adres.Gdata = den / ad1.Gdata;
		adres.Pdata = den / ad1.Pdata;
		adres.Adata = den / ad1.Adata;
		adres.user = den / ad1.user;
		return adres;
	}

	// division (*AziData)
	AziData& operator/=( const AziData& ad2 ) { 
		azi /= ad2.azi; Gdata /= ad2.Gdata; 
		Pdata /= ad2.Pdata; Adata /= ad2.Adata;
		user /= ad2.user;
		return *this; 
	}
	friend AziData operator/( const AziData& ad1, const AziData& ad2 ) {
		AziData adres = ad1;
		adres /= ad2;
		return adres;
	}

	// square root
	AziData sqrt() const { return AziData( std::sqrt(azi), std::sqrt(Gdata), std::sqrt(Pdata), std::sqrt(Adata), std::sqrt(user) ); }
	friend AziData sqrt( const AziData& ad1 ) { return ad1.sqrt();	}

	// fabs
	AziData fabs() const { return AziData( std::fabs(azi), std::fabs(Gdata), std::fabs(Pdata), std::fabs(Adata), std::fabs(user) ); }
	friend AziData fabs( const AziData& ad1 ) { return ad1.fabs();	}

	/* stream */
   friend std::ostream& operator<< ( std::ostream& o, const AziData& ad ) {
      o<<ad.azi<<"  "<<ad.Gdata<<"  "<<ad.Pdata<<"  "<<ad.Adata;
      return o;
   }

   static constexpr float NaN = -12345.;

};


/* data structure for measurements and predictions on a single station
	longitude, latitude, disance, azimuth, 
	measurements, GCP traveltimes, source terms, for groupT, phaseT, and amplitude */
struct StaData : public AziData {
   //bool valid = false;
   float lon = NaN, lat = NaN;
   float dis = NaN; //, azi = NaN;
	//float Gdata = NaN, Pdata = NaN, Adata = NaN;
	float Gsource = NaN, Psource = NaN, Asource = NaN;
	float Gpath = NaN, Ppath = NaN;

   StaData() : AziData() {}

/*
   StaData( float lonin, float latin, float disin, float Gdatain,
       float Pdatain, float aziin, float misGin, float misPin, float ampin )
      : AziData( aziin, misGin, misPin, ampin )
      , lon(lonin), lat(latin), dis(disin), Gdata(Gdatain), Pdata(Pdatain) {}
*/

	StaData( const float azi, const float lon, const float lat, const float Gdata, const float Pdata, const float Adata,
				const float init = NaN, const float userin = NaN )
		: lon(lon), lat(lat), Gpath(init), Gsource(init), Ppath(init), Psource(init), Asource(init),
		  AziData(azi, Gdata, Pdata, Adata, init) {
		if( userin != NaN ) user = userin;
	}

   StaData( const char *input, float ph_shift = 0. )
      : StaData() {
      sscanf(input, "%f %f %f %f %f", &lon, &lat, &Gdata, &Pdata, &Adata);
      Pdata += ph_shift;
   }

	// check if (all attributes are) within a given range defined by two other AziDatas
	virtual bool isWithin( const AziData& lbound, const AziData& ubound ) const {
		AziData ad;
		ToMisfit( ad );
		//ad.Pdata -= per * floor(ad.Pdata/per+0.5);
		return ad.isWithin( lbound, ubound );
	}

   bool isComplete() const { return ( lon!=NaN && lat!=NaN && Gdata!=NaN && Pdata!=NaN && Adata!=NaN ); }
   bool GisComplete() const { return ( Gdata!=NaN && Gpath!=NaN && Gsource!=NaN ); }
   bool PisComplete() const { return ( Pdata!=NaN && Ppath!=NaN && Psource!=NaN ); }
   bool AisComplete() const { return ( Adata!=NaN && Asource!=NaN ); }

	AziData ToAziData() const {
		return AziData( azi, Gdata, Pdata, Adata );
	}

	bool ToMisfit( AziData& ad ) const {
		if( ! ( GisComplete() && PisComplete() && AisComplete() ) )
			return false;
		ad.azi = azi;
		ad.Gdata = Gdata-Gpath-Gsource;
		ad.Pdata = Pdata-Ppath-Psource;
		ad.Adata = Adata-Asource;
		ad.user = Asource;
		return true;
	}

/*
   void streamTo ( std::ofstream& fout, int nformat = 1 ) const {
      if( nformat == 1 ) {
         fout<<azi<<"   "<<stdG<<" "<<misG<<"   "
             <<stdP<<" "<<misP<<"   "<<A<<" "<<Apred<<"\n";
      } else if( nformat == 2 ) {
         fout<<lon<<" "<<lat<<"   "<<Gdata-stdG<<" "<<stdG<<"   "
             <<Pdata-stdP<<" "<<stdP<<"   "<<Apred<<" "<<A<<"   "
             <<dis<<"   "<<azi<<"\n";
      }
   }
*/

   friend bool operator<( const StaData& ad1, const StaData& ad2 ) {
      return (ad1.azi < ad2.azi);
   }

	/*
   friend std::ostream& operator<< ( std::ostream& o, const StaData& sd ) {
      o<<"( "<<sd.lon<<", "<<sd.lat<<" ): dis = "<<sd.dis
       <<" Gdata = "<<sd.Gdata<<" Pdata = "<<sd.Pdata<<"    "<<sd.azi;
      return o;
   } */
   friend std::ostream& operator<< ( std::ostream& o, const StaData& sd ) {
      o<<sd.lon<<" "<<sd.lat<<"  "<<sd.dis<<" "<<sd.azi<<"  "
		 <<sd.Gdata<<" "<<sd.Gpath<<" "<<sd.Gsource<<"  "
		 <<sd.Pdata<<" "<<sd.Ppath<<" "<<sd.Psource<<"  "
		 <<sd.Adata<<" "<<sd.Asource;
      return o;
   }

protected:
   static constexpr float NaN = -12345.;

};


class ADAdder {
public:
	ADAdder( const bool useG, const bool useP, const bool useA ) {
		if(!useG && !useP && !useA) {
			fptr = &N;
		} else if( useG && !useP && !useA ) {
			fptr = &G;
		} else if( !useG && useP && !useA ) {
			fptr = &P;
		} else if( !useG && !useP && useA ) {
			fptr = &A;
		} else if( useG && useP && !useA ) {
			fptr = &GP;
		} else if( useG && !useP && useA ) {
			fptr = &GA;
		} else if( !useG && useP && useA ) {
			fptr = &PA;
		} else if( useG && useP && useA ) {
			fptr = &GPA;
		}
	}

	inline static float N( const AziData& ad ) { return 0; }
	inline static float G( const AziData& ad ) { return ad.Gdata; }
	inline static float P( const AziData& ad ) { return ad.Pdata; }
	inline static float A( const AziData& ad ) { return ad.Adata; }
	inline static float GP( const AziData& ad ) { return ad.Gdata+ad.Pdata; }
	inline static float GA( const AziData& ad ) { return ad.Gdata+ad.Adata; }
	inline static float PA( const AziData& ad ) { return ad.Pdata+ad.Adata; }
	inline static float GPA( const AziData& ad ) { return ad.Gdata+ad.Pdata+ad.Adata; }

	inline float operator()( const AziData& ad ) {
		return fptr(ad);
	}

private:
	float (*fptr)( const AziData& ad ) = nullptr;
};

#endif
