#ifndef SDCONTAINER_H
#define SDCONTAINER_H

#include "DataTypes.h"
#include "StackTrace.h"
#include "Map.h"
#include "RadPattern.h"
//#include "ModelInfo.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdexcept>

#ifndef FuncName
#define FuncName __FUNCTION__
#endif

/* exceptions */
namespace ErrorSC {

   class Base : public std::runtime_error {
   public:
      Base(const std::string message)
         : runtime_error(message) {
			PrintStacktrace();
      }
   };

   class BadFile : public Base {
   public:
      BadFile(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Cannot access file ("+info+").") {}
   };

   class BadParam : public Base {
   public:
      BadParam(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Bad parameters ("+info+").") {}
   };

   class InternalException : public Base {
   public:
      InternalException(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Internal exception ("+info+").") {}
   };
}

/* -------------------- data type enum -------------------- */
enum Dtype { Undefined=0, B, R, L }; // type of data to be used


/* ----- station data container ----- */
class SDContainer {
public:
	const float per = NaN;
	//const char type = 'N';
	const Dtype type = Undefined;

	/* con/destructors */

	// waveform: No vel maps needed. Misfits will be pushed back (to Pdata and Adata) on the fly. Path and source terms will be set to zero
	SDContainer( const float perin, const Dtype typein )
		: per(perin), type(typein) {
		if( type!=R && type!=L ) {
			std::string str( "unknown data type: " + std::to_string(type) ); //str.push_back(type);
			throw ErrorSC::BadParam( FuncName, str );
		}
	}

	// 3D model: Path velocities (and thus Tpaths) are predicted on the fly from vel maps
	SDContainer( const float perin, const Dtype typein, const std::string fmeasure, 
					 const std::string fmapG, const std::string fmapP, const std::string fsta="" ) 
		: per(perin), oop(1./perin), type(typein) {
		if( type!=R && type!=L ) {
			std::string str( "unknown data type: " + std::to_string(type) ); //str.push_back(type);
			throw ErrorSC::BadParam( FuncName, str );
		}
		LoadMeasurements( fmeasure, fsta );
		LoadMaps( fmapG, fmapP );
	}

	// 1D model: Path velocities are fixed at the input velocities
	SDContainer( const float perin, const Dtype typein, const std::string fmeasure, 
					 const float velG, const float velP, const std::string fsta="" ) 
		: per(perin), oop(1./perin), type(typein), 
		  _velG(velG), _velP(velP) {
		if( type!=R && type!=L ) {
			std::string str( "unknown data type: " + std::to_string(type) ); //str.push_back(type);
			throw ErrorSC::BadParam( FuncName, str );
		}
		LoadMeasurements( fmeasure, fsta );
	}
	/*
	SDContainer( float perin, const std::vector<StaData>& datain )
		: per(perin), dataV(datain) {}

	SDContainer( float perin, std::vector<StaData>&& datain )
		: per(perin), dataV( std::move(datain) ) {}
	*/

	void Sort() { std::sort(dataV.begin(), dataV.end()); }

	void push_back( const StaData& sd ) { dataV.push_back(sd); }

	StaData& back() { return dataV.back(); }

	std::size_t size() const { return dataV.size(); }

	/* dump into an AziData vector */
	//void ToAziVector( std::vector<AziData>& adV );
	void ToMisfitV( std::vector<AziData>& adV, const float Tmin ) const;

	/* compute azimuth and distance to a given center location */
	void UpdateAziDis( const float srclon, const float srclat );
	// predict/store traveltimes from VelMaps into Gpath&Ppath
	bool UpdatePathPred( const float srclon, const float srclat, const float srct0 );
	// predict/store source terms into Gpath&Ppath
	void UpdateSourcePred( const RadPattern& );
	// scale source amplitudes to match the observations
	void ComputeAmpRatios( std::vector<float>& ampratioV ) const {
		for( const auto& sd : dataV )
         if( sd.Adata!=NaN && sd.Asource!=NaN ) 
				ampratioV.push_back(sd.Adata/sd.Asource);
	}
	void AmplifySource( const float Afactor ) {
		for( auto& sd : dataV )
			if( sd.Asource != NaN ) sd.Asource *= Afactor;
	}
	// correct 2 pi for phase T misfits
	void Correct2PI() {
		for( auto& sd : dataV )
			sd.Pdata -= per * floor( (sd.Pdata-sd.Ppath-sd.Psource) *oop+0.5);
	}

	/* compute bin average */
	void BinAverage( std::vector<AziData>& adVmean, std::vector<AziData>& adVvar, const bool c2pi = true, const bool isFTAN = true );
	void BinAverage_ExcludeBad( std::vector<StaData>& sdVgood, const bool c2pi = true );

	/* IO */
	void LoadMeasurements( const std::string& fmeasure, const std::string fsta="" );
	void LoadMaps( const std::string& fmapG, const std::string& fmapP );
	void PrintAll( std::ostream& sout = std::cout ) {
		for( const auto& sd : dataV )
			sout<<sd<<"\n";
	}


protected:
   /* define class scope constants */
   static constexpr int MIN_BAZI_SIZE = 1;						/* allowed min number of measurements in each azimuth bin */
   static constexpr float BINSTEP = 20;							/* bin averaging step size */
   static constexpr float BINHWIDTH = 10;							/* bin averaging half width */

	static constexpr float exfactor = 2.5;							/* #sigma for excluding bad data */

   static constexpr float Min_Perc = 0.8;							/* allowed min fraction of path length in the vel map */
   static constexpr float Lfactor = 2.5;							/* define lamda = per * Lfactor for PathAvg */

   static constexpr float DISMIN = 0.;								/* allowed min */
   static constexpr float DISMAX = 9999.;							/* and max event-station distance for location searching */

   static constexpr float stdGest = 4.5;							/* an estimation of GroupT (sec), */
   static constexpr float stdPest = 1.5;							/* PhaseT (sec), */
   static constexpr float stdPHest = 0.5;							/* Spectrum Phase Shift (radians) */
   static constexpr float stdAest = 0.35;							/* and Amplitude (as fraction of the amplitude!) std-dev */
   static constexpr float varRGmin = 4.0, varLGmin = 4.0;		/* the lowerbound of GroupT (sec), was 0.8 */
   static constexpr float varRPmin = 0.5, varLPmin = 0.5;		/* PhaseT (sec), was 0.3 */
   static constexpr float varRPHmin = 0.05, varLPHmin = 0.05;	/* Specturm Phase Shift (radians), */
   static constexpr float varRAmin = 0.04, varLAmin = 0.04;		/* and Amplitude (as fraction of the amplitude square!) std-devs, was 0.02 */



   static constexpr int pio4_R = 0;									/* for initial phase test. normaly be 0. */
   static constexpr int pio4_L = 0;									/* is added to the FTAN coefficient piover4 */

	static constexpr float NaN = AziData::NaN;

private:
	const float oop = NaN;

	//ModelInfo model;	// bad idea! not necessary and adds inter-module complexity
	float lon = NaN, lat = NaN, t0 = NaN, stk = NaN, rak = NaN, dip = NaN, dep = NaN;

	Map mapG, mapP;
	const float _velG = NaN, _velP = NaN;
	std::vector<StaData> dataV;

	void HandleBadBins(std::vector<AziData>& adVmean, std::vector<AziData>& adVstd, const AziData adest ) const;
	// compute variance by propagating the given variance into the data
	void ComputeVariance( std::vector<AziData>& adVmean, std::vector<AziData>& adVstd, const AziData ad_varin ) const;
	// ( old: pull up any std-devs that are smaller than defined minimum )
	//void WaterLevel( std::vector<AziData>& adVmean, std::vector<AziData>& adVstd, const AziData ad_stdmin ) const;

};

#endif
