#ifndef EQKANALYZER_H
#define EQKANALYZER_H

#include "DataTypes.h"
#include "SDContainer.h"
#include "FileName.h"
#include <vector>
#include <map>

/* -------------------- data structures -------------------- */
enum Dtype { Undefined=0, B, R, L }; // type of data to be used

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
class EQKAnalyzer {
public:
	// model parameters
	ModelInfo _model;

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

   /* output misfit-v.s.-focal_corrections
    * should always be called after UpdateAziDis() and UpdateFocalCorr() */
   enum OutType { FIT, MAP };
   //void Output( bool excludeBad = false, const OutType otype = FIT );
   /* compute and output misfits, separately, for group, phase, and amplitudes */
   //void ComputeMisfitsAll();

	std::vector<float> perRlst() const, perLlst() const;
	void InitRadPattern();
	void chiSquare( const ModelInfo& minfo, float& chiS, float& wSum, int& N, bool isInit = false ) const;

private:
	// RadPattern objects for predicting source terms
	RadPattern _rpR, _rpL;

	//std::vector< std::vector<FileName> > fRlist, fLlist;
	/* store measurements/predictions of each station with a StaData,
	 * all StaDatas at a single period is handeled by a SDContainer */
	std::vector<SDContainer> _dataR, _dataL;

	/* ---------- input parameters ---------- */
	// search area of epicenter
	float clon, clat, ct0 = 0.;
	float Rs, Rlon, Rlat, Rt = 20.;	// Rlon: Rs(km) in longitude(deg)
	// data type
	char datatype_name;
	Dtype datatype;
	bool _useG = true, _useP = true, _useA = true;
	// data weightings
   float weightR_Loc, weightL_Loc;  // weighting between Rayleigh and Love data for Location search
   float weightR_Foc, weightL_Foc;  // weighting between Rayleigh and Love data for Focal search
	// Input files. Three files at each period for 1) measurements, 2) group vel map, and 3) phase vel map
	std::map<float, std::array<FileName, 3> > fRlist, fLlist;
	// Input eigen-function and phase-velocity files
	FileName fReigname, fRphvname;
   FileName fLeigname, fLphvname;
	// output files
	FileName outname_misF;           // filename for output focal misfit
   FileName outname_misL;           // filename for output location misfit
   FileName outname_misAll;         // filename for output all separated misfits (group, phase, amplitude)
   FileName outname_pos;            // filename for output posterior distribution
	std::map<float, FileName> outlist_RF, outlist_LF; // filename for output focal_fit
   std::map<float, FileName> outlist_RP, outlist_LP; // filename for output travel time predictions

private:
	// private functions
	inline void NormalizeWeights( Dtype& datatype, float& wR, float& wL);
	bool InitEpic();
	bool MKDir(const char *dirname);
	void MKDirFor( const std::string& path );

};

#endif
