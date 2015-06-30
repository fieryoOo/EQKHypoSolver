#ifndef SYNGENERATOR_H
#define SYNGENERATOR_H

#include "SacRec.h"
#include "ModelInfo.h"
#include <string>

class fstring : public std::string {
public:
	fstring() : std::string() {}
	fstring( const char* strin ) : std::string(strin) {}
	fstring( const std::string& strin ) : std::string(strin) {}
	fstring( const fstring& fstr2 ) : std::string(fstr2), pfstring(nullptr) {}
	~fstring() { clearfstring(); }

	char* f_str() const { return f_str( size()+1 ); }
	char* f_str( const size_t fsize ) const {
		if( fsize < size()+1 )
			throw std::runtime_error("Error(fstring::f_str): input size too small");
		clearfstring();
		pfstring = new char[fsize];
		sprintf(pfstring, "%s", c_str());
		std::fill(pfstring+size(), pfstring+fsize, ' ');
		return pfstring;
	}

	void assignf( const char* strin, const size_t strsize ) {
		int i=strsize-1;
		for(; i>=0; i--)
			if( strin[i] != ' ' ) break;
		this->assign( strin, i+1 );
	}

private:
	mutable char* pfstring = nullptr;
	void clearfstring() const { if(pfstring!=nullptr) delete [] pfstring; pfstring=nullptr; }
};


class SynGeneratorData {
protected:
	// internal data (used by fortran subroutines)
	// initial info
	fstring name_fmodel, name_feigen;
	bool key_compr = false;
	char sigR = '-', sigL = '-', its = 'C', modestr[2];
	//char fstr_fmodel[256];
	int im = 6, iq = 1, nd = 1000, nper; //npoints;
   float fix_vel = 2.8;
	float elat, elon, elatc, elonc, aM, tm[6], permin, permax, vmax = 100000; //dt;
	
	// station info
	int nsta = 0;
	char sta[2000][8], net[2000][8];
	float lat[2000], latc[2000], lon[2000];
   // surf_disp data
   float freq[2000], cr[2000], ur[2000], wvr[2000], ul[2000], cl[2000], wvl[2000];
   float ampr[2000], ampl[2000], ratio[2000], qR[2000], qL[2000], I0[2000];
   float v[2000][3], dvdz[2000][3];
   //std::unique_ptr<float[]> pv( new float[2000*3]() ); float *v = pv.get();
   //std::unique_ptr<float[]> pdvdz( new float[2000*3]() ); float *dvdz = pdvdz.get();

	// tracer data (managed by unique_ptr)
	//std::unique_ptr<float[]> pcor;
	bool traced = false;
};

class SynGenerator : SynGeneratorData {
public:
	SynGenerator() {}
	SynGenerator( const fstring& name_fmodel, const fstring& name_fphvel, const fstring& name_feigen, const char wavetype, const int mode ) {
		Initialize( name_fmodel, name_fphvel, name_feigen, wavetype, mode );
	}
	SynGenerator( const SynGenerator& sg2 ) 
		: SynGeneratorData(sg2), minfo(sg2.minfo) {
		if( traced ) {
			size_t ncor = 2000*2*500;
			pcor.reset( new float[ncor]() );
			if( ! pcor )
				throw std::runtime_error("new failed for pcor!");
			std::copy(sg2.pcor.get(), sg2.pcor.get()+ncor, pcor.get());
		}
	}


	void Initialize( const fstring& name_fmodel, const fstring& name_fphvel, const fstring& name_feigen, const char wavetype, const int mode );

	// station list
	void LoadSta( const std::string name_fsta );
	void ClearSta() { nsta = 0; }
	void PushbackSta( const SacRec& sac );

	// event info
	void SetEvent( const ModelInfo mi );

	// produce synthetic as sac file
	bool ComputeSyn( const std::string& staname, const float slon, const float slat, int npts, float delta, 
						  float f1, float f2, float f3, float f4, SacRec& sacZ, SacRec& sacN, SacRec& sacE );

	//bool Synthetic( const float lon, const float lat, const std::string& chname,
	//					 const float f1, const float f2, const float f3, const float f4, SacRec& sac );

protected:
	static constexpr float GEO = 0.993277, pio180 = M_PI / 180.;

private:
	ModelInfo minfo;

	// tracer data managed by unique_ptr
	std::unique_ptr<float[]> pcor;

	// trace all event-station GC paths
	void TraceAll();

	//void GetParams( const std::string name_fparam );
	void ReadPerRange( const std::string& name_fphvel, const int mode );

	static void WriteSACHeader( SAC_HD& shd, const int npts, const float dt, const float elat, const float elon,
										 const std::string& sta, const float slat, const float slon, 
										 const std::string& chn, const std::string& net );
};


#endif
