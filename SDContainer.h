#ifndef SDCONTAINER_H
#define SDCONTAINER_H

#include "DataTypes.h"
#include "StackTrace.h"
#include "Map.h"
#include "ModelInfo.h"
#include "RadPattern.h"
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

   class InternalException : public Base {
   public:
      InternalException(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Internal exception ("+info+").") {}
   };
}

class SDContainer {
public:
	const float per = NaN;

	/* con/destructors */
	SDContainer( const float perin, const std::string fmeasure, 
					 const std::string fmapG, const std::string fmapP ) 
		: per(perin), oop(1./perin) {
		LoadMeasurements( fmeasure );
		LoadMaps( fmapG, fmapP );
	}

	/*
	SDContainer( float perin, const std::vector<StaData>& datain )
		: per(perin), dataV(datain) {}

	SDContainer( float perin, std::vector<StaData>&& datain )
		: per(perin), dataV( std::move(datain) ) {}
	*/

	std::size_t size() const { return dataV.size(); }

	/* dump into an AziData vector */
	//void ToAziVector( std::vector<AziData>& adV );
	void ToMisfitV( std::vector<AziData>& adV ) const;

	/* compute azimuth and distance to a given center location */
	void UpdateAziDis( const float srclon, const float srclat );
	// predict/store traveltimes from VelMaps into Gpath&Ppath
	void UpdatePathPred( const float srclon, const float srclat, const float srct0 );
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
	void BinAverage( std::vector<AziData>& adVmean, std::vector<AziData>& adVstd );
	void BinAverage_ExcludeBad( std::vector<StaData>& sdVgood );

	/* IO */
	void LoadMeasurements( const std::string& fmeasure );
	void LoadMaps( const std::string& fmapG, const std::string& fmapP );
	void PrintAll( std::ostream& sout = std::cout ) {
		for( const auto& sd : dataV )
			sout<<sd<<"\n";
	}

	// for debug!
	friend int main( int argc, char* argv[] );

protected:
   /* define class scope constants */
   static constexpr int MIN_BAZI_SIZE = 1;      /* allowed min number of measurements in each azimuth bin */
   static constexpr float BINSTEP = 20;         /* bin averaging step size */
   static constexpr float BINHWIDTH = 10;       /* bin averaging half width */

	static constexpr float exfactor = 2.5;			/* #sigma for excluding bad data */

   static constexpr float Min_Perc = 0.95;      /* allowed min fraction of path length in the vel map */
   static constexpr float Lfactor = 2.;         /* define lamda = per * Lfactor for PathAvg */

   static constexpr float DISMIN = 0.;          /* allowed min */
   static constexpr float DISMAX = 9999.;       /* and max event-station distance for location searching */

   static constexpr float stdGest = 4.5;        /* an estimation of GroupT, */
   static constexpr float stdPest = 1.5;        /* PhaseT, */
   static constexpr float stdAest = 0.42;       /* and Amplitude (as fraction of the amplitude!) std-dev */

   static constexpr float stdGmin = 3.0;        /* the lowerbound of GroupT, was 0.8 */
   static constexpr float stdPmin = 1.0;			/* PhaseT, was 0.3 */
   static constexpr float stdAmin = 0.05;			/* and Amplitude (as fraction of the amplitude square!) std-devs, was 0.02 */

   static constexpr float Pthreshold = 0.005;   /* the threshold for probability in searching for parameter
                                                   sensitivity prior to the Monte Carlo search, */
   static constexpr float Sfactor = 0.1;        /* and the step half-length for the search as a fraction
                                                   of (ub-lb) decided by Pthreshold */

   static constexpr int pio4_R = 0;             /* for initial phase test. normaly be 0. */
   static constexpr int pio4_L = 0;             /* is added to the FTAN coefficient piover4 */

	static constexpr float NaN = AziData::NaN;

private:
	const float oop = NaN;
	//float srclon = NaN, srclat = NaN;
	ModelInfo model;
	Map mapG, mapP;
	std::vector<StaData> dataV;

	void HandleBadBins(std::vector<AziData>& adVmean, std::vector<AziData>& adVstd, const AziData adest ) const;
	void WaterLevel( std::vector<AziData>& adVmean, std::vector<AziData>& adVstd, const AziData ad_stdmin ) const;

};

#endif
