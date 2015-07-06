#include "SynGenerator.h"
#include "SacRec.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

/* FORTRAN entrance */

#ifdef __cplusplus
extern"C" {
#endif
//	void read_params_( char name_feigen[255], char name_fparam[255], int* mode, char name_fmodel[256], float* elat, float* elon, float* aM,
//							 float* dt, float tm[6], int* im, int* iq, char* its, int* nd, int* npoints, int* nper, float* tmin, float* tmax, float* vmax,
//							 int* nsta, char sta[2000][8], float lat[2000], float latc[2000], float lon[2000], char net[2000][8] );

	void angles2tensor_(float* stk, float* dip, float* rak, float tm[6]);

	void atracer_( char name_fmodel[256], float* elon, float* elat, int* ncor, bool* applyQ, int* nsta,
						float latc[2000], float lon[2000], float* cor );

//	void surfread_( char name_feigen[255], char* sigR, char* sigL, char modestr[2], int* nper, int* nd, float* depth, float freq[2000],
	void surfread_( char *feig_buff, int *eig_len, char* sigR, char* sigL, char modestr[2], int* nper, int* nd, float* depth, float freq[2000],
						float cr[2000], float ur[2000], float wvr[2000], float cl[2000], float ul[2000], float wvl[2000], float v[2000][3], 
						float dvdz[2000][3], float ampr[2000], float ampl[2000], float ratio[2000], float qR[2000], float qL[2000], float I0[2000] );

	void cal_synsac_( int* ista, char* its, char* sigR, char* sigL, float* cor, float* f1, float* f2, float* f3, float* f4,
							float* vmax, float* fix_vel, int* iq, int* npoints, float freq[2000], float* dt, int* nper,
							bool* key_compr, float* elatc, float* elonc, float qR[2000], float qL[2000], int* im, float* aM, float tm[6],
							float ampl[2000], float cr[2000], float ul[2000], float ur[2000], float wvl[2000], float wvr[2000],
							float v[2000][3], float dvdz[2000][3], float ratio[2000], float I0[2000], float* slatc, float* slon,
							float* sigz, float* sign, float* sige );

#ifdef __cplusplus
}
#endif


void SynGenerator::Initialize( const fstring& name_fmodel_in, const fstring& name_fphvel, const fstring& name_feigen_in, const char wavetype, int mode ) {
	name_fmodel = name_fmodel_in;
	name_feigen = name_feigen_in;
	// input params
	// wave type
	//const char wavetype = 'R';
	switch( wavetype ) {
		case 'R': sigR = '+'; break;
		case 'L': sigL = '+'; break;
		default: throw std::runtime_error("Error(SynGenerator::SynGenerator): invalid wave type input!");
	}
	// mode number
	mode += 1;	// the fortran code takes mode=1 as the fundamental mode
	sprintf(modestr, "%d", mode);
	if( modestr[1] == '\0' ) modestr[1] = ' ';
	/*
	if(argc==7) {
		key_compr = true;
		fix_vel = atof(argv[6]);
	}
	*/

	//std::cout<<"Wavetype = "<<wavetype<<" mode# = "<<mode<<std::endl;

/*
	// call read_params
	read_params_( name_feigen.f_str(255), name_fparam.f_str(255), &mode, fstr_fmodel, &elat, &elon, &aM, &dt, 
					  tm, &im, &iq, &its, &nd, &npoints, &nt, &tmin, &tmax, &vmax,
					  &nsta, sta, lat, latc, lon, net );
	elatc = atan(0.993277*tan(elat*M_PI/180.)); elonc = elon * M_PI/180.;
*/

	// update permin, permax, and nt by fphvel
	ReadPerRange( name_fphvel, mode );

	// read feigen into memory
	std::ifstream fin( name_feigen );
	if( ! fin ) throw std::runtime_error("Error(Initialize): IO failed on "+name_feigen);
	fin.seekg(0, std::ios::end); feig_len = fin.tellg();
	peig.reset( new char[feig_len] );
	fin.seekg(0, std::ios::beg);
	fin.read(peig.get(), feig_len);
}

void SynGenerator::ReadPerRange( const std::string& name_fphvel, const int mode ) {
	// open file
	std::ifstream fin(name_fphvel);
	if( ! fin )
		throw std::runtime_error("IO failed on " + name_fphvel);
	// read first two lines
	float perlast = -123.; permin = -123.;
	std::string line;
	if( std::getline(fin, line) ) sscanf(line.c_str(), "%f", &permin);
	if( std::getline(fin, line) ) sscanf(line.c_str(), "%f", &perlast);
	if( permin==-123. || perlast==-123. )
		throw std::runtime_error("invalid phavel file: " + name_fphvel);
	const float dper = perlast - permin;
	// read till last
	fin.clear(); fin.seekg(0);
	int nmode = 0;
	bool validlastline = false;
	nper = 2;
	while( std::getline(fin, line) ) {
		bool validline = (sscanf(line.c_str(), "%f", &permax) == 1);
		if( validline && !validlastline ) { // new mode starts
			if( nmode == mode ) break;	// break if the last mode is wanted
			nmode++;
			nper = 0;
		}
		if( validline && validlastline ) {
			if( dper != permax-perlast )
				throw std::runtime_error("invalid phavel file: " + name_fphvel);
		}
		if( validline ) nper++;
		validlastline = validline;
		perlast = permax;
	}
	if( nmode != mode ) // mode not found!
		throw std::runtime_error("mode number not found in " + name_fphvel);
	if( nper > 1998 )	// memory will explode
		throw std::runtime_error("too many periods in " + name_fphvel);
	//std::cerr<<"nper = "<<nper<<" permin = "<<permin<<" permax = "<<permax<<" for mode "<<mode<<std::endl;

}


void SynGenerator::LoadSta( const std::string name_fsta ) {
	// update station info
	std::ifstream fin(name_fsta);
	if( ! fin )
		throw std::runtime_error("IO failed on " + name_fsta);
	nsta = 0;
	for( std::string line; std::getline(fin,line); ) {
		if( sscanf(line.c_str(), "%s %s %f %f", net[nsta], sta[nsta], &(lat[nsta]), &(lon[nsta])) != 4 )
			continue;
		latc[nsta] = atan(GEO*tan(pio180*lat[nsta]))/pio180;
		nsta++;
	}

	// needs re-trace
	traced = false;
}

void SynGenerator::PushbackSta( const SacRec& sac ) {
	strcpy(sta[nsta], sac.stname().c_str());
	strcpy(net[nsta], sac.ntname().c_str());
	lon[nsta] = sac.shd.stlo;
	lat[nsta] = sac.shd.stla;
	latc[nsta] = atan(GEO*tan(pio180*lat[nsta]))/pio180;
	//std::cerr<<sta[nsta]<<" "<<net[nsta]<<" "<<lon[nsta]<<" "<<lat[nsta]<<" "<<latc[nsta]<<std::endl;
	nsta++;
	traced = false;
}

void SynGenerator::SetEvent( const ModelInfo mi ) {
	//const float GEO = 0.993277, pio180 = M_PI / 180.;
	// update (internal) model info
	minfo = mi;
	// update event location
	elat = mi.lat; elon = mi.lon;
	elatc = atan(GEO*tan(elat*pio180)); elonc = elon * pio180;
	// update focal mechanism
	aM = mi.M0;
	angles2tensor_(&minfo.stk, &minfo.dip, &minfo.rak, tm);

	// call surfread with the new depth
	//#pragma omp critical
	surfread_( peig.get(), &feig_len, &sigR, &sigL, modestr, &nper, &nd, &(minfo.dep), freq, cr, ur, wvr,
				  cl, ul, wvl, v, dvdz, ampr, ampl, ratio, qR, qL, I0 );
	//surfread_( name_feigen.f_str(255), &sigR, &sigL, modestr, &nper, &nd, &(minfo.dep), freq, cr, ur, wvr,

	// needs re-trace
	traced = false;
}

void SynGenerator::TraceAll() {
	// call atracer
	int ncor;
	bool applyQ = true;
   //std::unique_ptr<float[]> pcor( new float[2000*2*500]() ); float *cor = pcor.get();
	pcor.reset( new float[2000*2*500]() );
	if( ! pcor )
		throw std::runtime_error("new failed for pcor!");
	//#pragma omp critical
	atracer_( name_fmodel.f_str(256), &elat, &elon, &ncor, &applyQ, &nsta, latc, lon, pcor.get() );
	//atracer_( fstr_fmodel, &elat, &elon, &ncor, &applyQ, &nsta, latc, lon, pcor.get() );
	// all traced
	traced = true;
}

bool SynGenerator::ComputeSyn( const std::string& staname, const float slon, const float slat, int npts, float delta,
										 float f1, float f2, float f3, float f4, SacRec& sacz, SacRec& sacn, SacRec& sace ) {
	// trace all GCPs
	if( ! traced ) TraceAll();

	/*/ calc base size
	int nbase = 2; n2pow = 1;
	while( n2pow<13 && nbase<npoints ) { n2pow++; nbase <<= 1; }
	std::cerr<<"n2pow = "<<n2pow<<"  nbase = "<<nbase<<"  npoints = "<<npoints<<"\n";
	nf = nbase / 2;
	df = 1. / (dt*nbase);// f0 = 0.; */

	// period range
   //if(TMAX.gt.(nt-1)*dper+per1)TMAX=PERMAX-2.*dper	// these are what misha had in the fortran version
   //if(TMIN.lt.per1)TMIN=PERMIN+2.*dper					// cannot understand the logic...
	if( 1./f1>permax || 1./f4<permin )
		throw std::runtime_error("filter corner freq out of range!");

	// search for the requested station in the list
	int ista = 0;
	for(; ista<nsta; ista++) {
		if( staname == std::string(sta[ista]) &&
			 slon == lon[ista] && slat == lat[ista] ) break;
	}
	if( ista == nsta ) return false;

	// station found. Produce synthetics
	std::string netname = net[ista];// staname = sta[ista];
	//std::cerr<<ista<<" STA = "<<staname<<" NET = "<<netname<<" STLAT = "<<lat[ista]<<" STLON = "<<lon[ista]<<"\n";
	// prepare sac headers
	//SacRec sacz, sacn, sace;
	std::string chnprefix = delta>=1 ? "LH" : "BH";
	auto init_sac = [&]( const char comp, SacRec& sac ) {
		std::string chn = chnprefix + comp;
		WriteSACHeader( sac.shd, npts, delta, elat, elon, staname, lat[ista], lon[ista], chn, netname );
		sac.ResizeSig();
	};
	init_sac( 'Z', sacz ); init_sac( 'N', sacn ); init_sac( 'E', sace );
	//sacz.MutateAs(sacz); sacn.MutateAs(sacz); sace.MutateAs(sacz);
	int ista_f = ista+1;
	cal_synsac_( &ista_f, &its, &sigR, &sigL, pcor.get(), &f1, &f2, &f3, &f4, &vmax, &fix_vel, &iq,
			&npts, freq, &delta, &nper, &key_compr, &elatc, &elonc, qR, qL,
			&im, &aM, tm, ampl, cr, ul, ur, wvl, wvr, v, dvdz, ratio, I0,
			&(latc[ista]), &(lon[ista]), sacz.sig.get(), sacn.sig.get(), sace.sig.get() );
	// flip (misha's coordinate is upside down)
	sacz.Mul(-1.); sacn.Mul(-1.); sace.Mul(-1.);
	// shift b times to match the data (origin-time = b-time-of-real-data + t0)
	sacz.shd.b += minfo.t0;	sacn.shd.b += minfo.t0;	sace.shd.b += minfo.t0;

	return true;
}


void SynGenerator::WriteSACHeader( SAC_HD& shd, const int npts, const float delta, const float elat, const float elon,
														const std::string& sta, const float slat, const float slon, 
														const std::string& chn, const std::string& net ) {
	shd.npts = npts;
	shd.delta = delta;
	shd.b = 0.; 
	shd.e = shd.delta*(shd.npts-1);
	shd.evdp = 0.0;
	shd.evel = 0.0;
	shd.evla = elat;
	shd.evlo = elon;
	strcpy(shd.kstnm,sta.c_str());
	shd.stdp = 0.;
	shd.stel = 0.;
	shd.stla = slat;
	shd.stlo = slon;
	shd.nzyear = 2000;
	shd.nzjday = 1;
	shd.nzhour = shd.nzmin = shd.nzsec = shd.nzmsec = 0;
	strcpy(shd.kcmpnm, chn.c_str());
	strcpy(shd.knetwk, net.c_str());
}



/*
bool SynGenerator::Synthetic( const float lon, const float lat, const std::string& chname, 
										const float f1, const float f2, const float f3, const float f4, SacRec& sac ) {
	std::string saclist("sacS.list");
	bool found = false;
	std::ifstream fin(saclist);
	if( ! fin )
		throw std::runtime_error("cannot access file "+saclist);
	for( std::string line; std::getline(fin, line); ) {
		std::stringstream ss(line); ss >> line;
		sac.Load( line );
		if( chname==sac.chname() && lon==sac.shd.stlo && lat==sac.shd.stla ) {
			found = true;
			break;
		}
	}
	if( ! found ) {
		sac.clear();
	} else {
		sac.Mul(-1.);	// is the synthetic upside down???
		sac.Filter(f1,f2,f3,f4);
		sac.shd.b -= minfo.t0;
		sac.Resample();
	}
	return found;
}
*/

