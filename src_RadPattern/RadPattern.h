#ifndef RADPATTERN_H
#define RADPATTERN_H

#include "MyOMP.h"
#include "EigenRec.h"
#include <cmath>
#include <memory>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>


/* ---------- exceptions ---------- */
//#define FuncName __PRETTY_FUNCTION__
#define FuncName __FUNCTION__
namespace ErrorRP {
   class BadFile : public std::runtime_error {
   public:
      BadFile(const std::string funcname, const std::string info = "")
	 : runtime_error("Error("+funcname+"): Cannot access file ("+info+").") {}
   };

   class BadParam : public std::runtime_error {
   public:
      BadParam(const std::string funcname, const std::string info = "")
        : runtime_error("Error("+funcname+"): Bad parameters ("+info+").") {}
   };

   class BadAzi : public std::runtime_error {
   public:
      BadAzi(const std::string funcname, const std::string info = "")
        : runtime_error("Error("+funcname+"): Unexpected azimuths ("+info+").") {}
   };

   class BadBuff : public std::runtime_error {
   public:
      BadBuff(const std::string funcname, const std::string info = "")
        : runtime_error("Error("+funcname+"): Internal buffer modified ("+info+"). The RadPattern class is not yet thread-safe!") {}
   };

   class HeaderMismatch : public std::runtime_error {
   public:
      HeaderMismatch(const std::string funcname, const std::string info = "")
        : runtime_error("Error("+funcname+"): HeaderMismatch ("+info+").") {}
   };

};


typedef float ftype;
typedef std::map<float, std::array<float,3>> MA3;

class RadPattern {
public:
	char type; 
	float dep, M0;
	std::array<float, 6> MT;
	//float stk, dip, rak;

	RadPattern( const char type = 'N', const std::string& feigname = "" );

   void SetModel( const char type, const std::string& feigname );

	friend std::ostream& operator <<(std::ostream &o, const RadPattern &rp) {
		o<<rp.type<<" "<<rp.grtM.size()<<"   "
		 <<rp.dep<<" "<<rp.M0<<" "<<rp.MT[0]<<" "<<rp.MT[1]<<" "<<rp.MT[2]<<" "<<rp.MT[3]<<" "<<rp.MT[4]<<" "<<rp.MT[5];
		 //<<rp.stk<<" "<<rp.dip<<" "<<rp.rak<<" "<<rp.dep<<" "<<rp.M0;
		return o;
	}

   /* Predict Rayleigh/Love wave radiation patterns */
   bool Predict( const ftype strike, const ftype dip, const ftype rake,
					  const ftype depth, const ftype M0, const std::vector<float>& perlst ) {
		MT = MomentTensor(strike, dip, rake); Predict(MT, depth, M0, perlst);
	}
   bool Predict( const std::array<ftype, 6>& MT, const ftype depth, const ftype M0, const std::vector<float>& perlst );

	/* get the amp norm term at per */
	std::array<float, 2> cAmp( const float per ) const;

	/* prediction at one single azimuth. return false if the given azimuth is invalidated due to small amplitude */
	// M0 = scalar seismic momentum
	// dis = distance; alpha = attenuation coeff
	// J = mode energy integration (from eigen);
	// U = local group velocity
	bool GetPred( const float per, const float azi,	float& grt, float& pht, float& amp, const float mul = 1. ) const;
	bool GetPred( const float per, const float azi,	float& grt, float& pht, float& amp,
					  const float dis, const float alpha, const float recCAmp = NaN ) const;

	void OutputPreds( const std::string& fname, const float norm_dis = -1, const float Q = -1 ) const;

	RadPattern& operator-=(const RadPattern &rp2);
	friend RadPattern operator-(const RadPattern &rp1, const RadPattern &rp2) {
		auto rp3 = rp1; rp3 -= rp2; return rp3;
	}
	RadPattern& operator*=(const float&);

	// re-scale (the amplitude of) *this/rp1s to have the minimum amp chiSquare misfit to rp2
	// sigmaV is a map (by period) of sigmaAs (fraction 0-1)
	void NormBy( const RadPattern &rp2, const MA3 &sigmasM );
	void NormBy( const RadPattern &rp2, const std::map<float,float> &sigmaM );
	friend void NormAmps( RadPattern &rp1R, RadPattern &rp1L, const RadPattern &rp2R, const RadPattern &rp2L,
								 const MA3 &sigmasMR, const MA3 &sigmasML );

	// compute the chi square misfits of G,P,and A between *this and rp2 ( returns {chiSG, chiSP, chiSA, n})
	// *this or rp1R&rp1L is/are modified when normalizedA==true
	// sigmasV is a map (by period) of sigma triples: sigmaG (sec), sigmaP (sec), sigmaA (fraction 0-1)
	std::array<float,4> chiSquare(const RadPattern &rp2, const MA3 &sigmasM, const bool normalizeA=true);
	friend std::array<float,4> chiSquare(RadPattern &rp1R, RadPattern &rp1L, const RadPattern &rp2R, const RadPattern &rp2L,
													 const MA3 &sigmasMR, const MA3 &sigmasML, const bool normalizeA=true);

public:
	static constexpr float NaN = -12345.;
	//static constexpr float oofourpi = 0.25e-13 / M_PI;	// unit convertion = 1.0e-13
	static constexpr int nazi = 181;
	static constexpr int dazi = 2;
	static constexpr int InvalidateHwidth = 3;	// half width in iazi of the focal pred invalidating window
	static constexpr float AmpValidPerc = 0.05;  // focal predictions with A < Aaverage*AmpValidPerc are invalidated

private:
   //struct Rimpl; std::unique_ptr<Rimpl> pimplR;
	EigenRec er;

	// sample azimuths
	std::vector<float> aziV;
	// I0 (mod energy integral) keyed by period
	//std::map< float, float > I0M;
	// campM[per][0]: source amp norm term : campM[per] = 1.e-15/( (phvel*grvel*I0) * sqrt(8 * pi) )
	// campM[per][1]: angular wave number
	std::map< float, std::array<float,2> > campM;
	// group, phase, amplitudes keyed by period
	std::map< float, std::vector<float> > grtM, phtM, ampM;

	
	std::array<float, 6> MomentTensor( float stk, float dip, float rak, const float M0=1. ) const;

	void ShiftCopy( std::vector<float>& Vout, const float* arrayin, const int nazi ) const;

	void NormCoefs( const RadPattern &rp2, const std::map<float,float> &sigmaM, float &a, float &b );

	//void LoadEig();
};

#endif
