#ifndef SACREC_H
#define SACREC_H

#include "mysac64.h"
#include "MyOMP.h"
//#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <limits>
#include <stdexcept>
#include <cmath>

#ifndef FuncName
#define FuncName __FUNCTION__
#endif
/*
// a workaround for nullptr when compiled by old gcc compiler
const class {				// this is a const object...
public:
   template<class T>			// convertible to any type
   operator T*() const { return 0 }	// of null non-member pointer...
   template<class C, class T>		// or any type of null
   operator T C::*() const { return 0; }// member pointer...
private:
  void operator&() const;		// whose address can't be taken
} nullptr = {};				// and whose name is nullptr
*/


namespace ErrorSR {

   class Base : public std::runtime_error {
   public:
		Base( const std::string funcname, const std::string message )
			: runtime_error(funcname + ": " + message) {
				//PrintStacktrace();
			}
	};

   class BadFile : public Base {
   public:
      BadFile(const std::string funcname, const std::string info = "")
         : Base(funcname, "Invalid or non-accessable file ("+info+").") {}
   };

   class BadParam : public Base {
   public:
      BadParam(const std::string funcname, const std::string info = "")
         : Base(funcname, "Bad parameters ("+info+").") {}
   };

   class EmptySig : public Base {
   public:
      EmptySig(const std::string funcname, const std::string info = "")
         : Base(funcname, "No sac signal loaded in the memory ("+info+").") {}
   };

   class InsufData : public Base {
   public:
      InsufData(const std::string funcname, const std::string info = "")
         : Base(funcname, "Insufficient data points ("+info+").") {}
   };

   class HeaderMismatch : public Base {
   public:
      HeaderMismatch(const std::string funcname, const std::string info = "")
         : Base(funcname, "Header mismatch ("+info+").") {}
   };

	
   class UndefMethod : public Base {
   public:
		UndefMethod(const std::string funcname, const std::string info = "")
         : Base(funcname, "Undefined method ("+info+").") {}
   };

   class ExternalError : public Base {
   public:
		ExternalError(const std::string funcname, const std::string info = "")
         : Base(funcname, "External error ("+info+").") {}
   };

   class MemError : public Base {
   public:
		MemError(const std::string funcname, const std::string info = "")
         : Base(funcname, "Memory error ("+info+").") {}
   };
};


class SacRec {
public:
   std::string fname;			// input file name
   SAC_HD shd;				// sac header
   std::unique_ptr<float[]> sig;	// pointer to the signal
   //std::auto_ptr<float> sig;
public:
   /* ------------------------------ con/destructors and operators ------------------------------ */
   /* constructors */
	SacRec( std::ostream& reportin = std::cerr );	// default 1
   SacRec( const std::string& fnamein, std::ostream& reportin = std::cerr );	// default 2
	SacRec( const size_t npts, std::ostream& reportin = std::cerr );	// default 3
   SacRec( const SacRec& recin );		// copy
   SacRec( SacRec&& recin );			// move
   /* operators */
   SacRec &operator= ( const SacRec& recin );	// assignment
   SacRec &operator= ( SacRec&& recin );	// move
   /* destructor */
   ~SacRec(); 

	/* check signal/header validation */
	void updateDeps();
	bool isValid() {
		updateDeps();
		return (shd.delta>0 && shd.npts>0 && shd.depmin==shd.depmin && shd.depmax==shd.depmax);
	}
	bool isZero() {
		float* sigsac = sig.get();
		for(int i=0; i<shd.npts; i++)	if( sigsac[i] != 0. ) return false;
		return true;
	}

	/* assign header and allocate memory */
	void MutateAs ( const SacRec& recin );

	/* allocate memory for signal */
	void ResizeSig() { ResizeSig(shd.npts); }
	void ResizeSig( const size_t npts ) {
		if( npts <= 0 )
			throw ErrorSR::BadParam( FuncName, "negative npts!");
		shd.npts = npts; shd.e = shd.b + shd.delta*npts;
		sig.reset(new float[npts]);
		if( ! sig )
			throw ErrorSR::MemError( FuncName, "new failed!");
	}

   /* ------------------------------ sac file read/write ------------------------------ */
   /* load sac header from file 'fname' */
   void LoadHD( const std::string& fnamein ) { fname = fnamein; LoadHD(); }
   void LoadHD();
   /* read sac header+signal from file 'fname', memory is allocated on heap */
   void Load( const std::string& fnamein ) { fname = fnamein; Load(); }
   void Load();
	/* clear sac and release memory */
	void clear() { sig.reset(); shd = sac_null; fname.clear(); }
   /* write to file '*fname' */
   void WriteHD( const std::string& fname );
   void Write( const std::string& fname );
	/* load a txt file */
	void LoadTXT( const std::string& fname );
	/* dump signal to stdout/txt */
	void Dump( const std::string fname = "" );
	/* dump header to stdout/txt */
	void DumpHD( const std::string fname = "" );
	/* print a single header field */
	void PrintHD( const std::string field, std::ostream &o = std::cout ) const;
   /* change a single header filed */
	void ChHdr(const std::string& field, const std::string& value){
		std::istringstream sin( field + " " + value ); sin >> shd;
	}

	float Dis() const;
	float Dis();
	float Azi() const;
	float Azi();

	const std::string ntname() const {
		std::stringstream ss(shd.knetwk);
		std::string ntname; ss >> ntname;
		return ntname;
	}
	const std::string evname() const {
		std::stringstream ss(shd.kevnm);
		std::string evname; ss >> evname;
		return evname;
	}
	const std::string stname() const {
		std::stringstream ss(shd.kstnm);
		std::string stname; ss >> stname;
		return stname;
	}
	const std::string chname() const {
		std::stringstream ss(shd.kcmpnm);
		std::string chname; ss >> chname;
		return chname;
	}

   /* ------------------------------ header/signal information ------------------------------ */
	inline size_t Index( const float time ) const;
	inline float Time( const size_t index ) const { return shd.b + index*shd.delta; }
   /* compute the absolute time in sec relative to 1900.01.00 */
   double AbsTime ();
   /* update/reformat header time if shd.nzmsec is modified and is out of the range [0,1000) */
   void UpdateTime();
   /* search for min&max signal positions and amplitudes */
	void MinMax (int& imin, int& imax) const { MinMax( imin, imax, shd.b, shd.e ); }
	void MinMax (int& imin, int& imax, float tbegin, float tend) const;
   void MinMax ( float tbegin, float tend, float& tmin, float& min, float& tmax, float& max ) const;
   /* compute the root-mean-square average in a given window */
   float RMSAvg ( float tbegin, float tend ) const { return RMSAvg( tbegin, tend, 1); }
   float RMSAvg ( float tbegin, float tend, int step ) const;
	bool Mean ( float& mean ) const { return Mean(shd.b, shd.e, mean); }
	bool Mean ( float tbegin, float tend, float& mean ) const { return Mean(tbegin, tend, 1, mean); }
	bool Mean ( float tbegin, float tend, int step, float& mean ) const;
	bool MeanStd ( float& mean, float& std ) const { return MeanStd(shd.b, shd.e, mean, std); }
	bool MeanStd ( float tbegin, float tend, float& mean, float& std ) const { return MeanStd(tbegin, tend, 1, mean, std); }
	bool MeanStd ( float tbegin, float tend, int step, float& mean, float& std ) const;
	// compute accurate time/amplitude of the peak (fit with a parabola)
	float Tpeak() const { float t, a; Peak(t, a); return t; }
	float Tpeak( const float tbegin, const float tend ) const { float t, a; Peak(t, a, tbegin, tend); return t; }
	float Apeak() const { float t, a; Peak(t, a); return a; }
	float Apeak( const float tbegin, const float tend ) const { float t, a; Peak(t, a, tbegin, tend); return a; }
	void Peak(float& tpeak, float& apeak) const { Peak(tpeak, apeak, shd.b, shd.e); }
	void Peak(float& tpeak, float& apeak, const float tbegin, const float tend) const;
	float Sig( float time ) const;	// compute accurate sig value at a given time (fit with a parabola)

   /* ------------------------------ single-sac operations ------------------------------ */
	template<class Functor>	void Transform(const Functor& func, size_t ib=0, int ie=NaN) {
		if( !sig ) throw ErrorSR::EmptySig(FuncName);

		if( ib < 0 ) ib = 0;
		if( ie==NaN || ie>shd.npts ) ie = shd.npts;
		float* sigsac = sig.get();
		for(int i=ib; i<ie; i++)	func(sigsac[i]);
	}

	template<class Functor>	void Transform2(const Functor& func, size_t ib=0, int ie=NaN) {
		if( !sig ) throw ErrorSR::EmptySig(FuncName);

		if( ib < 0 ) ib = 0;
		if( ie==NaN || ie>shd.npts ) ie = shd.npts;
		float* sigsac = sig.get();
		for(int i=ib; i<ie; i++)	func(Time(i), sigsac[i]);
	}

   void Mul( const float mul );
	void Addf( const SacRec& sac2 );
	void Subf( const SacRec& sac2 );
	void Divf( const SacRec& sac2 );

	void Reverse() { 
		SacRec sac2;
		Reverse(sac2);
		*this = std::move(sac2);
	}
	void Reverse( SacRec& sac2 );

	float SNR( const float tsignall, const float tsignalh, const float tnoisel, const float tnoiseh ) const;

	/* performs integration in the time domain using the trapezoidal rule */
	void IntegrateT() { IntegrateT(*this); }
	void IntegrateT( SacRec& sac_out ) const;

	/* performs integration in the frequency domain (omega arithmetic) */
	void Integrate() { Integrate(*this); }
	void Integrate( SacRec& sac_out ) const;

	/* performs differentiation in the frequency domain (omega arithmetic) */
	void Differentiate() { Differentiate(*this); }
	void Differentiate( SacRec& sac_out ) const;

	void PullUpTo( const SacRec& sac2 );
   void ToAm() { ToAm(*this);	}
   void ToAm( SacRec& sac_am ) const {
		SacRec sac_ph;
		ToAmPh( sac_am, sac_ph );
	}
	void FFT( SacRec& sac_re, SacRec& sac_im, const int nfout = 0 ) const;	// in series when sig is large
	void FFT_p( SacRec& sac_re, SacRec& sac_im, const int nfout = 0 ) const;	// always parallel
	void ToAmPh( SacRec& sac_am, SacRec& sac_ph, const int nfout = 0 ) const;	// in series when sig is large
	void ToAmPh_p( SacRec& sac_am, SacRec& sac_ph, const int nfout = 0 ) const;	// always parallel
	void FromAmPh( SacRec& sac_am, SacRec& sac_ph, const short outtype = 0 );		// in series when sig is large
	void FromAmPh_p( SacRec& sac_am, SacRec& sac_ph, const short outtype = 0 );	//	always parallel
	/* filters */
	void LowpassCOSFilt( double f3, double f4 ) { LowpassCOSFilt(f3, f4, *this); }
	void LowpassCOSFilt( double f3, double f4, SacRec& srout ) { Filter(-1., -1., f3, f4, 0, srout); }
	void HighpassCOSFilt( double f1, double f2 ) { HighpassCOSFilt(f1, f2, *this); }
	void HighpassCOSFilt( double f1, double f2, SacRec& srout ) { Filter(f1, f2, -1., -1., 1, srout); }
	void BandpassCOSFilt( double f1, double f2, double f3, double f4 ) { BandpassCOSFilt(f1, f2, f3, f4, *this); }
	void BandpassCOSFilt( double f1, double f2, double f3, double f4, SacRec& srout ) { Filter(f1, f2, f3, f4, 2, srout); }
	void LowpassBTWFilt( double fc, int n ) { LowpassBTWFilt(fc, n, *this); }
	void LowpassBTWFilt( double fc, int n, SacRec& srout ) { Filter(-1., -1., fc, n, 3, srout); }
	void HighpassBTWFilt( double fc, int n ) { HighpassBTWFilt(fc, n, *this); }
	void HighpassBTWFilt( double fc, int n, SacRec& srout ) { Filter(fc, n, -1., -1., 4, srout); }
	void BandpassBTWFilt( double fcL, double fcR, int n, bool zeroPhase = false ) {
		BandpassBTWFilt(fcL, fcR, n, *this, zeroPhase); 
	}
	void BandpassBTWFilt( double fcL, double fcR, int n, SacRec& srout, bool zeroPhase = false ) { 
		Filter(-1., fcL, fcR, n, 5, srout, zeroPhase); 
		/*
		if( zeroPhase ) {
			srout.Reverse();
			srout.Filter(-1., fcL, fcR, n, 5);
			srout.Reverse();
		} 
		*/
	}
	void GaussianFilt( double fc, double fhlen ) { GaussianFilt(fc, fhlen, *this); }
	void GaussianFilt( double fc, double fhlen, SacRec& srout ) { Filter(-1., fc, fhlen, -1., 6, srout); }
   /* method that performs different types of filters:
	 * type = 0: Lowpass cosine -f3~f4_
	 * type = 1: highpass cosine _f1~f2-
	 * type = 2: bandpass cosine _f1~f2-f3~f4_
	 * type = 3: lowpass butterworth -fc=f3~n=f4_
	 * type = 4: highpass butterworth _fc=f1~n=f2-
	 * type = 5: bandpass butterworth _fcL=f1~nL=f2-fcR=f3~nR=f4_
	 * type = 6: gaussian _fc=f2~fhlen=f3_ */
   void Filter ( double f1, double f2, double f3, double f4, const int type, bool zeroPhase = false ) {	// in-place (in series when sig is large)
		Filter(f1, f2, f3, f4, type, *this, zeroPhase); 
	}
   void Filter ( double f1, double f2, double f3, double f4, const int type, SacRec& srout, bool zeroPhase = false );	// out-of-place (in series when sig is large)
   void Filter_p ( double f1, double f2, double f3, double f4, const int type, SacRec& srout, bool zeroPhase = false );	// always parallel
	/* tapers */
	void cosTaperL( const float fl, const float fh );
	void cosTaperR( const float fl, const float fh );
	void gauTaper( const float fc, const float fh );
   /* remove mean and trend */
   void RTrend();
   /* remove response and apply filter 
		type: 0=displacement, 1=velocity, 2=acceleration */
   void RmRESP( const std::string& fresp, float perl, float perh, const int type = 1 ) {
		std::string evrexe;
		RmRESP( fresp, perl, perh, evrexe, type );
	}
   void RmRESP( const std::string& fresp, float perl, float perh, const std::string& evrexe, const int type = 1 );
   /* resample (with anti-aliasing filter) the signal to given sps */
   //void Resample( bool fitParabola = true ) { Resample( floor(1.0/shd.delta+0.5), fitParabola ); }
   void Resample( int sps = -1, bool fitParabola = true );
	/* smoothing ( running average ) */
	void Smooth( float timehlen, SacRec& sacout ) const;
	void Hilbert() { Hilbert(*this); }
	void Hilbert( SacRec& sacout );
	void Envelope() { Envelope(*this); }
	void Envelope( SacRec& sacout );

   void cut( float tb, float te ) { cut(tb, te, *this); }
   void cut( float tb, float te, SacRec& );
   /* ------------------------------ inter-sac operations ------------------------------ */
   /* merge a second sacrec to the current */
   void Merge( SacRec sacrec2 ) {
      merge( sacrec2 );
      arrange();
   }
   void merge( SacRec sacrec2 );
   int arrange( const char *recname = nullptr );

	/* ---------- compute the correlation coefficient with an input SacRec ---------- */
	float Correlation( const SacRec& sac2 ) const {
		return Correlation( sac2, shd.b, shd.e );
	}
	float Correlation( const SacRec& sac2, const float tb, const float te ) const;

	/* Cross-Correlate with another sac record
		ctype=0: Cross-Correlate (default) 
		ctype=1: deconvolve (sac.am/sac2.am)
		ctype=2: deconvolve (sac2.am/sac.am) */
	void CrossCorrelate( SacRec& sac2 ) { CrossCorrelate(sac2, *this); }
	void CrossCorrelate( SacRec& sac2, SacRec& sacout, int ctype = 0 );

	/* ------------------------------- cut by event ---------------------------------- */
	void ZoomToEvent( const std::string etime, float evlon, float evlat, float tb, float tlen, std::string ename = "" );
	void ZoomToEvent( const SAC_HD& eshd, float evlon, float evlat, float tb, float tlen, std::string ename );

	/* ------------------------------- temporal normalizations ------------------------------- */
	void OneBit();
	void RunAvg( float timehlen, float Eperl, float Eperh );

	/* ------------------------------- memory consumed ------------------------------- */
	float MemConsumed() const;
	void AlwaysParallel() { maxnpts4parallel = std::numeric_limits<int>::max(); }
	// to run the fftw, 16 times the original npts is required ( in&out complex double array with size doubled for specturm ). 20 is used to be safe
	void SetMaxMemForParallel( float MemInMb ) { maxnpts4parallel = (MemInMb * 1024. * 1024. - 1000.) / (4. * 20.); }

	// define string output
   friend std::ostream& operator<< ( std::ostream& o, const SacRec& sac ) {
		o << sac.fname << "  " << sac.evname() << " " << sac.shd.evlo << " " << sac.shd.evla
		  << "  " << sac.stname() << " " << sac.shd.stlo << " " << sac.shd.stla;
		return o;
	}

	static constexpr float NaN = -12345.;

protected:
	int maxnpts4parallel = 1e6;

private:
   /* impl pointer */
   struct SRimpl;
   std::unique_ptr<SRimpl> pimpl;
	/* reporting stream */
	std::ostream* report = &(std::cerr);

};

/*
class AMPH {
private:
   class APimpl;
   std::unique_ptr<APimpl> pimpl;
   std::unique_ptr<float[]> am, ph;
};
*/

#endif
