#ifndef EQKANALYZER_H
#define EQKANALYZER_H

#include "ModelInfo.h"
#include "SynGenerator.h"
#include "RadPattern.h"
#include "SDContainer.h"
#include "Searcher.h"
#include "FileName.h"
#include "SacRec.h"
#include <vector>
#include <array>
#include <map>

//class ModelInfo;class SDContainer;

/* -------------------- data type -------------------- */
//enum Dtype { Undefined=0, B, R, L }; // type of data to be used

/* -------------------- exceptions -------------------- */
namespace WarningEA {
   class Base {
   public:
      Base(const std::string message) {
         if( nWarnOther < nWarnMax ) {
            std::cerr<<message<<std::endl;
            nWarnOther++;
         } else if( nWarnOther == nWarnMax ) {
            std::cerr<<"### Too many warnings. Further warning outputs truncated! ###"<<std::endl;
            nWarnOther++;
         }
      }
   private:
      static const int nWarnMax = 2e5;
      static int nWarnOther;
   };

   class MoveExistFile : public Base {
   public:
      MoveExistFile(const std::string funcname, const std::string info = "")
         : Base("Warning("+funcname+"): Moving existing file ("+info+").") {}
   };

   class BadFile : public Base {
   public:
      BadFile(const std::string funcname, const std::string info = "")
         : Base("Warning("+funcname+"): Unable to access file ("+info+").") {}
   };

   class BadParam : public Base {
   public:
      BadParam(const std::string funcname, const std::string info = "")
         : Base("Warning("+funcname+"): Bad parameters ("+info+").") {}
   };

   class Other : public Base {
   public:
      Other(const std::string funcname, const std::string info = "")
         : Base("Warning("+funcname+"): "+info) {}
   };
};

namespace ErrorEA {

   class Base : public std::runtime_error {
   public:
      Base(const std::string message)
         : runtime_error(message) {
      }
		~Base() { PrintStacktrace(); }
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

   class SizeMismatch : public Base {
   public:
      SizeMismatch(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Incompatible sizes ("+info+").") {}
   };

   class BadAzi : public Base {
   public:
      BadAzi(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Bad azi triplet ("+info+").") {}
   };

   class EmptyData : public Base {
   public:
      EmptyData(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Empty data input ("+info+").") {}
   };

   class InsufData : public Base {
   public:
      InsufData(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Insufficient data points ("+info+").") {}
   };

   class BadPred : public Base {
   public:
      BadPred(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Invalid focal prediction ("+info+").") {}
   };

   class InternalException : public Base {
   public:
      InternalException(const std::string funcname, const std::string info = "")
         : Base("Error("+funcname+"): Internal exception ("+info+").") {}
   };
};


/* -------------------- the analyzer class -------------------- */
class EQKAnalyzer : public Searcher::IDataHandler<ModelInfo> {
public:
	// model parameters
	//ModelInfo _model;

   /* con/destructors and operators */
   EQKAnalyzer();
   EQKAnalyzer( const std::string, bool MoveExistF = true );
   //EQKAnalyzer( const EQKAnalyzer& );
   //EQKAnalyzer( EQKAnalyzer&& );
   //EQKAnalyzer& operator= ( const EQKAnalyzer& );
   //EQKAnalyzer& operator= ( EQKAnalyzer&& );
   //~EQKAnalyzer();

   /* prepare database */
   void LoadParams( const FileName&, const bool MoveExistF = true );
   int Set( const char*, const bool MoveExistF = true );
   void CheckParams();
   void LoadData();
	void SaveOldOutputs() const;

	inline std::vector<float> perRlst() const;
	inline std::vector<float> perLlst() const;
	inline std::vector<float> perlst(const Dtype) const;

	// initialize the Analyzer by pre- predicting radpatterns and updating pathpred for all SDContainer
	// version (1): non-const, modifies internal data
	void UpdatePredsM( const ModelInfo& mi ) { UpdatePredsM( mi, _dataR, _dataL ); }
	// version (2): const, modify external data 
	void UpdatePredsM( const ModelInfo& mi, std::vector<SDContainer> &dataR, std::vector<SDContainer> &dataL ) const;

	void UpdatePredsW( const ModelInfo& minfo, std::vector<SDContainer> &dataR, std::vector<SDContainer> &dataL ) const;
	void UpdatePredsW( const ModelInfo& minfo ) {
		UpdatePredsW( minfo, _dataR, _dataL );
	}

	void UpdatePreds( const ModelInfo& minfo ) {
		if( _usewaveform ) {
			UpdatePredsW( minfo );
		} else {
			UpdatePredsM( minfo );
		}
	}

/*
	void FilldataW( const ModelInfo& minfo ) {
		float chiS; int N;
		SDContainer dataR, dataL;
		chiSquareW( minfo, chiS, N, true, dataR, dataL );
		_dataR.clear(); _dataL.clear(); 
		if(dataR.type!=Undefined) _dataR.push_back( std::move(dataR) );
		if(dataL.type!=Undefined) _dataL.push_back( std::move(dataL) );
	}
*/

	// use Love group data only when isInit=true
	void SetInitSearch( bool isInit ) { _isInit = isInit; }

	// compute M0 for least chiS when true
	void SetCorrectM0( bool correct ) { _correctM0 = correct; }

	// chi-square misfits from measurements-predictions
	void chiSquareM( const ModelInfo &minfo, float& chiS, int& N ) const;

	// chi-square misfits from waveform data-synthetics
	void chiSquareW( const ModelInfo &minfo, float& chiS, int& N ) const;
/*
	void chiSquareW( ModelInfo minfo, float& chiS, int& N, bool filldata, SDContainer& dataR, SDContainer& dataL ) const;
	void chiSquareW( ModelInfo minfo, float& chiS, int& N ) const {
		SDContainer dataR, dataL;
		chiSquareW( minfo, chiS, N, false, dataR, dataL );
	}
*/

	// call the relevant chiSquare method based on the _usewaveform value
	inline void chiSquare( const ModelInfo& minfo, float& chiS, int& N ) const {
		if( _usewaveform ) {
			chiSquareW( minfo, chiS, N );
		} else {
			chiSquareM( minfo, chiS, N );
		}
	}

	void Energy( const ModelInfo& minfo, float& E, int& Ndata ) const {
		float chiS; 
		chiSquare( minfo, chiS, Ndata );
		E = chiS * _indep_factor; //chiS/(Ndata-8.);
		if( Ndata < NdataMin )
			throw ErrorEA::InsufData( FuncName, std::to_string(Ndata) + " < " + std::to_string(NdataMin) );
		//return isvalid;
	}

   /* output misfit-v.s.-focal_corrections
    * should always be called after UpdateAziDis() and UpdateFocalCorr() */
   //enum OutType { FIT, MAP };
   //void Output( bool excludeBad = false, const OutType otype = FIT );
   //void ComputeMisfitsAll();

	void EstimateSigmas( const ModelInfo& minfo );
	void OutputSigmas() const;
   // output G,P,A data and predictions to separated files for each period
	void OutputFits( const ModelInfo &minfo );
   // compute and output misfits, separately, for group, phase, and amplitudes
	void OutputMisfits( const ModelInfo &minfo );
	// output source predictions (continuously in azimuth, for group, phase, and amplitudes) into single file for R/L waves
	void OutputSourcePatterns( const ModelInfo& mi );
	// output real (processed) and synthetic waveforms when the waveform fitting method is used
	void OutputWaveforms( const ModelInfo& mi ) { OutputWaveforms( mi, outdir_sac );	}
	void OutputWaveforms( const ModelInfo& mi, const std::string& outdir );

public:
	/* ---------- input parameters that needs to be externally accessible ---------- */
	float _indep_factor = 1.;			// describes the correlations among data (0: 100% correlated, 1: 100% independent)
	FileName outname_misF;           // filename for output focal misfit
   FileName outname_misL;           // filename for output location misfit
   FileName outname_misAll;         // filename for output all separated misfits (group, phase, amplitude)
   FileName outname_pos;            // filename for output posterior distribution
	FileName outname_srcR;				// filename for output Rayl source patterns
	FileName outname_srcL;				// filename for output Love source patterns
	FileName outdir_sac;					// dirname for output real and synthetic sacs
	FileName outname_sigmas;			// filename for output sigmas (predefined or estimated)

protected:
	static const bool dynamicVar = false;			// when true, uncertainties are dynamically estimated in each bin for computing chiSquare 
	static const int NdataMin = 3;
	static constexpr float NaN = AziData::NaN;
	static const bool rotateSyn = false;			// false=R&T; true=N&E
   static constexpr float _SNRMIN = 18.;			// allowed min SNR for waveform fitting */
   static constexpr float _DISMIN = 0.;			// allowed min */
   static constexpr float _DISMAX = 2000.;		// and max event-station distance for location searching */

private:	// class
	struct SinglePeriodInfo {
		float sigmaG = -1, sigmaP = -1, sigmaA = -1;
		FileName fmeasure, fmapG, fmapP, fstalst;
	};

private: // variables
	const size_t nthd = omp_get_max_threads();
	// option 1. measurements
	// RadPattern objects for predicting source terms (a pair for each thread)
	mutable std::vector<RadPattern> _rpR{std::vector<RadPattern>(nthd)}, _rpL{std::vector<RadPattern>(nthd)};

	/* store measurements/predictions of each station with a StaData,
	 * all StaDatas at a single period is handeled by a SDContainer */
	std::vector<SDContainer> _dataR, _dataL;

	// option 2. waveform fitting data
	typedef std::array<SacRec, 3> SacRec3;
	std::vector<SacRec3> _sac3VR, _sac3VL;
	SynGenerator _synGR, _synGL;

	// initial location ( for computing dis when loading data )
	float initlon, initlat;

	/* ---------- input parameters ---------- */
	// search area of epicenter
	//float clon, clat, ct0 = 0.;
	//float Rs, Rlon, Rlat, Rt = 20.;	// Rlon: Rs(km) in longitude(deg)
	// data selection
   float SNRMIN = _SNRMIN;	// for waveform fitting method only 
	float DISMIN = _DISMIN, DISMAX = _DISMAX;		
	// data type
	char datatype_name;
	Dtype datatype;
	bool _useG = true, _useP = true, _useA = true;
	bool _usewaveform = false;
	bool _isInit = false;
	// data weightings (!!!not implemented, adjust varmins in SDContainer instead!!!)
   float weightR_Loc = 1., weightL_Loc = 1.;  // weighting between Rayleigh and Love data for Location search
   float weightR_Foc = 1., weightL_Foc = 1.;  // weighting between Rayleigh and Love data for Focal search
	// Input files. Three files at each period for 1) measurements, 2) group vel map, and 3) phase vel map
	//std::map<float, std::array<FileName, 4> > fRlist, fLlist;
	std::map<float, SinglePeriodInfo> spiRM, spiLM;
	// Input eigen-function and phase-velocity files
	FileName fReigname, fRphvname;
   FileName fLeigname, fLphvname;
	// Input saclists and model file for the waveform synthetic
	FileName fmodelR, fmodelL;
	FileName fsaclistR, fsaclistL;
	short sacRtype = NaN, sacLtype = NaN;	// expecting 0 for displacement or 1 for velocity
	float f1 = NaN, f2 = NaN, f3 = NaN, f4 = NaN;
	/* ---------- internal variables ---------- */
	//bool _source_updated = false;
	bool _correctM0 = false;
	//float _AfactorR = 1, _AfactorL = 1;
	// output files
	std::map<float, FileName> outlist_RF, outlist_LF; // filename for output focal_fit
   //std::map<float, FileName> outlist_RP, outlist_LP; // filename for output travel time predictions

private: // functions
	bool FilenameToVel( const FileName& fname, float& vel ) const;
	inline void NormalizeWeights( Dtype& datatype, float& wR, float& wL);
	inline float ShiftInto( float val, float lb, float ub, float T) const;
	inline float BoundInto( float val, float lb, float ub ) const;
	//bool InitEpic();
	bool MKDir(const char *dirname) const;
	void MKDirs( const std::string& path ) const { MKDirFor(path, true); }
	void MKDirFor( const std::string& path, const bool isdir = false ) const;

	//float Tpeak( const SacRec& sac ) const;
	SacRec ComputeSyn(const SacRec &sac, SynGenerator &synG) const;
	StaData WaveformMisfit( const SacRec3 &sac3, SynGenerator &synG ) const;
	float RescaleSourceAmps( std::vector<SDContainer>& dataR, std::vector<SDContainer>& dataL ) const;
};

#endif
