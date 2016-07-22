#include "EQKAnalyzer.h"
#include "SynGenerator.h"
#include "Parabola.h"
//#include "VectorOperations.h"
//#include "DataTypes.h"
#include <sstream>
#include <sys/stat.h>
#include <cstdlib>

int WarningEA::Base::nWarnOther = 0;

EQKAnalyzer::EQKAnalyzer() {}

EQKAnalyzer::EQKAnalyzer( const std::string fparam, bool MoveF ) {
	LoadParams( fparam, MoveF );
	CheckParams();
}

//EQKAnalyzer::~EQKAnalyzer() {}



/* -------------------- param/data preparations (loading) -------------------- */
bool EQKAnalyzer::MKDir(const char *dirname) const {
	//create dir if not exists
	//with read/write/search permissions for owner and group, and with read/search permissions for others if not already exists
	if( mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0 ) return true;
	switch(errno) {
		case EEXIST:
			return false;
		default:
			perror("### Error: MKDir failed"); //failed. prompt to continue
			exit(0);
	}
}

void EQKAnalyzer::MKDirFor( const std::string& path, bool isdir ) const {
	if( path.empty() ) return;
	std::stringstream spath(path);
	std::vector<std::string> dirs;
	for(std::string dir; std::getline(spath, dir, '/'); ) {
		dirs.push_back(dir);
	}
	int dsize = dirs.size();
	if( ! isdir ) dsize--;
	std::string pathcur;
	for( int idir=0; idir<dsize; idir++ ) {
		pathcur += dirs[idir];
		if( MKDir( pathcur.c_str() ) )
			WarningEA::Other(FuncName, "Making directory: " + pathcur );
		pathcur += "/";
	}
}

void EQKAnalyzer::LoadParams( const FileName& fname, const bool MoveF ) {
   fname.CheckAccess();
   std::ifstream fin(fname);
   if( ! fin ) throw ErrorEA::BadFile(FuncName, fname);
   int nparam = 0;
   for( std::string stmp; std::getline(fin, stmp); ) {
      int retval = Set( stmp.c_str(), MoveF );
      if( retval == -3 ) continue; // invalid input for this parameter
      else if( retval == -2 ) continue; // empty input
      else if( retval == -1 ) continue;// std::cerr<<"Warning(EQKSearcher::Load): Unknown parameter name: "<<stmp<<std::endl;
      else if( retval == 0 ) continue;// std::cerr<<"Warning(EQKSearcher::Load): Empty parameter field for "<<stmp<<std::endl;
      else nparam++;
   }
   fin.close();

   std::cout<<"### EQKAnalyzer::LoadParams: "<<nparam<<" succed loads from param file "<<fname<<". ###\n"
				<<"    disRange = "<<DISMIN<<" - "<<DISMAX<<" perRange = "<<1./f3<<" - "<<1./f2<<"\n"
				<<"    SNRmin = "<<SNRMIN<<"(for waveform-fitting) datatype = "<<datatype_name<<" indep_factor = "<<_indep_factor<<std::endl;

   // check clon/clat validation and initialize einfo
   //if( ! InitEpic() ) throw ErrorEA::BadParam(FuncName, "invalid Rs | clon | clat");

   // synchronize output vectors with input-measurement vectors
   //pimplES->ArrangeOutlist( pimplES->outlist_RF, pimplES->perRtmp_fit, pimplES->perRlist );
   //pimplES->ArrangeOutlist( pimplES->outlist_LF, pimplES->perLtmp_fit, pimplES->perLlist );
   //pimplES->ArrangeOutlist( pimplES->outlist_RP, pimplES->perRtmp_pred, pimplES->perRlist );
   //pimplES->ArrangeOutlist( pimplES->outlist_LP, pimplES->perLtmp_pred, pimplES->perLlist );

   // check output paths and make directorie(s) if necessary
   std::vector<FileName> Voutf { outname_misF, outname_misL, outname_misAll, outname_pos };
	for( auto fname : outlist_RF ) { Voutf.push_back( fname.second+"_sta" ); Voutf.push_back( fname.second+"_bin" ); }
	for( auto fname : outlist_LF ) { Voutf.push_back( fname.second+"_sta" ); Voutf.push_back( fname.second+"_bin" ); }
	//for( auto fname : outlist_RP ) Voutf.push_back( fname.second );
	//for( auto fname : outlist_LP ) Voutf.push_back( fname.second );
   for( const auto& outfname : Voutf ) MKDirFor( outfname );
}

int EQKAnalyzer::Set( const char *input, const bool MoveF ) {
	std::istringstream buff(input);
	std::string stmp;
	if( ! (buff>>stmp) ) return -2;
	bool succeed;
	if( stmp == "fRse" ) succeed = (bool)(buff >> fReigname);
	else if( stmp == "fRsp" ) succeed = (bool)(buff >> fRphvname);
	else if( stmp == "fLse" ) succeed = (bool)(buff >> fLeigname);
	else if( stmp == "fLsp" ) succeed = (bool)(buff >> fLphvname);
	else if( stmp == "weightR_Loc" ) succeed = (bool)(buff >> weightR_Loc);
	else if( stmp == "weightL_Loc" ) succeed = (bool)(buff >> weightL_Loc);
	else if( stmp == "weightR_Foc" ) succeed = (bool)(buff >> weightR_Foc);
	else if( stmp == "weightL_Foc" ) succeed = (bool)(buff >> weightL_Foc);
	else if( stmp == "noG" ) { succeed = true; _useG = false; }
	else if( stmp == "noP" ) { succeed = true; _useP = false; }
	else if( stmp == "noA" ) { succeed = true; _useA = false; }
	else if( stmp == "SNRMIN" ) succeed = (bool)(buff >> SNRMIN);
	else if( stmp == "DISMIN" ) succeed = (bool)(buff >> DISMIN);
	else if( stmp == "DISMAX" ) succeed = (bool)(buff >> DISMAX);
	else if( stmp == "lon" ) { 
		succeed = (bool)(buff >> initlon); 
		if( succeed && initlon<0.) initlon += 360.; 
	}
	else if( stmp == "lat" ) succeed = (bool)(buff >> initlat);
	else if( stmp == "fmodelR" ) succeed = (bool)(buff >> fmodelR);
	else if( stmp == "fmodelL" ) succeed = (bool)(buff >> fmodelL);
	else if( stmp == "fsaclistR" ) {
		succeed = (bool)(buff >> fsaclistR >> sacRtype);
		if( sacRtype!=0 && sacRtype!=1 )
			throw ErrorEA::BadParam( FuncName, "sacRtype: expecting either 0 (displacement) or 1 (velocity)");
	}
	else if( stmp == "fsaclistL" ) {
		succeed = (bool)(buff >> fsaclistL >> sacLtype);
		if( sacLtype!=0 && sacLtype!=1 )
			throw ErrorEA::BadParam( FuncName, "sacLtype: expecting either 0 (displacement) or 1 (velocity)");
	}
	else if( stmp == "permin" ) {
		float permin;
		succeed = (bool)(buff >> permin);
		f3 = 1./permin; f4 = 1.1/permin;
	}
	else if( stmp == "permax" ) {
		float permax;
		succeed = (bool)(buff >> permax);
		f1 = 0.9/permax; f2 = 1./permax;
	}
	else if( stmp == "indep" ) { 
		succeed = (bool)(buff >> _indep_factor);
		if( _indep_factor<=0. || _indep_factor>1. )
			throw ErrorEA::BadParam( FuncName, "indep factor out of range: "+std::to_string(_indep_factor) );
	}
	else if( stmp == "dflag" ) {
		succeed = (bool)(buff >> datatype_name);
		if( succeed ) {
			switch( datatype_name ) {
				case 'B': datatype = B;
							 break;
				case 'R': datatype = R;
							 break;
				case 'L': datatype = L;
							 break;
				default: return -3;
			}
		}
	}
	else if( stmp == "fRm" ) {
		FileName stmp1, stmp2, stmp3, stmp4;
		float per;
		succeed = (bool)(buff>>stmp1>>stmp2>>stmp3>>per);	// fmeasure, fmapG, fmapP, per
		buff >> stmp4;													// station list (optional)
		if( succeed ) {
			//if( spiRM.find(per) != spiRM.end() ) return -3;	// period already exists
			auto &spi = spiRM[per];
			spi.fmeasure = stmp1; spi.fmapG = stmp2; spi.fmapP = stmp3; spi.fstalst = stmp4;
		}
	}
	else if( stmp == "fLm" ) {
		FileName stmp1, stmp2, stmp3, stmp4;
		float per;
		succeed = (bool)(buff>>stmp1>>stmp2>>stmp3>>per);	// fmeasure, fmapG, fmapP, per
		buff >> stmp4;													// station list (optional)
		if( succeed ) {
			//if( fLlist.find(per) != fLlist.end() ) return -3;	// period already exists
			auto &spi = spiLM[per];
			spi.fmeasure = stmp1; spi.fmapG = stmp2; spi.fmapP = stmp3; spi.fstalst = stmp4;
		}
	}
	else if( stmp == "sigmaR" ) {
		float per, sG, sP, sA;
		succeed = (bool)(buff>>per>>sG>>sP>>sA);
		if( succeed ) {
			auto &spi = spiRM[per];
			spi.sigmaG = sG; spi.sigmaP = sP; spi.sigmaA = sA;
		}
	}
	else if( stmp == "sigmaL" ) {
		float per, sG, sP, sA;
		succeed = (bool)(buff>>per>>sG>>sP>>sA);
		if( succeed ) {
			auto &spi = spiLM[per];
			spi.sigmaG = sG; spi.sigmaP = sP; spi.sigmaA = sA;
		}
	}
	else if( stmp == "fmisL" ) {
		FileName& outname = outname_misL;
		succeed = (bool)(buff >> outname);
		if( succeed && MoveF ) outname.SaveOld();
		/*
		if( succeed && MoveF && access(outname.c_str(), F_OK) == 0 ) {
			FileName oldname = outname + "_old";
			WarningEA::MoveExistFile(FuncName, outname+" -> "+oldname);
			rename(outname.c_str(), oldname.c_str());
		}
		*/
	}
	else if( stmp == "fmisF" ) {
		FileName& outname = outname_misF;
		succeed = (bool)(buff >> outname);
		if( succeed && MoveF ) outname.SaveOld();
	}
	else if( stmp == "fmisAll" ) {
		FileName& outname = outname_misAll;
		succeed = (bool)(buff >> outname);
		if( succeed && MoveF ) outname.SaveOld();
	}
	else if( stmp == "fpos" ) {
		FileName& outname = outname_pos;
		succeed = (bool)(buff >> outname);
		if( succeed && MoveF ) outname.SaveOld();
	}
	else if( stmp == "fsrcR" ) {
		FileName& outname = outname_srcR;
		succeed = (bool)(buff >> outname);
		if( succeed && MoveF ) outname.SaveOld();
	}
	else if( stmp == "fsrcL" ) {
		FileName& outname = outname_srcL;
		succeed = (bool)(buff >> outname);
		if( succeed && MoveF ) outname.SaveOld();
	}
	else if( stmp == "dirsac" ) {
		FileName& dirname = outdir_sac;
		succeed = (bool)(buff >> dirname);
		if( succeed && MoveF ) dirname.SaveOld();
	}
	else if( stmp == "ffitR" ) {
		FileName outname;
		float per;
		succeed = (bool)(buff >> outname >> per);
		if( succeed ) {
			outlist_RF[per] = outname;
			if( MoveF ) {
				FileName outname_sta = outname + "_sta";
				FileName outname_bin = outname + "_bin";
				outname_sta.SaveOld();	// move _sta file
				outname_bin.SaveOld();	// move _bin file
			}
		}
	}
	else if( stmp == "ffitL" ) {
		FileName outname;
		float per;
		succeed = (bool)(buff >> outname >> per);
		if( succeed ) {
			outlist_LF[per] = outname;
			if( MoveF ) {
				FileName outname_sta = outname + "_sta";
				FileName outname_bin = outname + "_bin";
				outname_sta.SaveOld();	// move _sta file
				outname_bin.SaveOld();	// move _bin file
			}
		}
	}
	else if( stmp == "fsigmas" ) {
		FileName& outname = outname_sigmas;
		succeed = (bool)(buff >> outname);
		if( succeed && MoveF ) outname.SaveOld();
	}
/*
	else if( stmp == "fpredR" ) {
		FileName outname;
		float per;
		succeed = (bool)(buff >> outname >> per);
		if( succeed ) {
			outlist_RP[per] = outname;
			//if( MoveF ) outname.SaveOld();
		}
	}
	else if( stmp == "fpredL" ) {
		FileName outname;
		float per;
		succeed = (bool)(buff >> outname >> per);
		if( succeed ) {
			outlist_LP[per] = outname;
			//if( MoveF ) outname.SaveOld();
		}
	} */
	else return -1;
	return succeed;
}


void EQKAnalyzer::SaveOldOutputs() const {
	outname_misF.SaveOld();
	outname_misL.SaveOld();
	outname_misAll.SaveOld();
	outname_pos.SaveOld();
	outname_srcR.SaveOld();
	outname_srcL.SaveOld();
	outname_sigmas.SaveOld();
	outdir_sac.SaveOld();
	for( const auto &pfpair : outlist_RF ) {
		((FileName)(pfpair.second + "_sta")).SaveOld();
		((FileName)(pfpair.second + "_bin")).SaveOld();
	}
	for( const auto &pfpair : outlist_LF ) {
		((FileName)(pfpair.second + "_sta")).SaveOld();
		((FileName)(pfpair.second + "_bin")).SaveOld();
	}
}

/* -------------------- param/data preparations (checking and post-processing) -------------------- */
inline void EQKAnalyzer::NormalizeWeights( Dtype& datatype, float& wR, float& wL) {
	if( datatype == Undefined ) datatype = B;
	//throw ErrorFS::BadParam(FuncName, "Undefined datatype");
	if( datatype == R ) {
		wR = 1.; wL = 0.;
	} else if( datatype == L ) {
		wR = 0.; wL = 1.;
	} else {
		if( wR<=0. || wL<=0. )
			throw ErrorEA::BadParam(FuncName, "invalid weightR | weightL");
		float factor = 2. / (wR + wL);
		wR *= factor;
		wL *= factor;
	}
}

void EQKAnalyzer::CheckParams() {
	// check initial focal info
	/*
	const float NaN = ModelInfo::NaN;
	if( _model.strike==NaN|| _model.dip==NaN||
		 _model.rake==NaN|| _model.depth==NaN) {
      throw ErrorEA::BadParam(FuncName, "invalid/incomplete initial focal info");
   }
	*/

   // check model inputs
   if( access(fReigname.c_str(), F_OK) == -1 ) throw ErrorEA::BadFile(FuncName, fReigname);
   if( access(fRphvname.c_str(), F_OK) == -1 ) throw ErrorEA::BadFile(FuncName, fRphvname);
   if( access(fLeigname.c_str(), F_OK) == -1 ) throw ErrorEA::BadFile(FuncName, fLeigname);
   if( access(fLphvname.c_str(), F_OK) == -1 ) throw ErrorEA::BadFile(FuncName, fLphvname);

	// normalize data weightings
	NormalizeWeights( datatype, weightR_Loc, weightL_Loc );
	NormalizeWeights( datatype, weightR_Foc, weightL_Foc );

	// check excluded data
   if( !_useG && !_useP && !_useA )
      throw ErrorEA::BadParam(FuncName, "All data excluded! (noG and noP)");
   if( ! _useG )
      WarningEA::Other(FuncName, "GroupT data excluded!" );
   if( ! _useP )
      WarningEA::Other(FuncName, "PhaseT data excluded!" );
   if( ! _useA )
      WarningEA::Other(FuncName, "Amplitude data excluded!" );

}


/* -------------------- Check/Load all data file requested by fparam -------------------- */
bool EQKAnalyzer::FilenameToVel( const FileName& fname, float& vel ) const {
	try {
		fname.CheckAccess();
		// succeed, is a file
		return false;
	} catch(...) {
		vel = strtof( fname.c_str(), NULL );
		if( vel<0.2 || vel>10. )
			throw ErrorEA::BadParam(FuncName, "input "+fname+" is neither a valid file nor a valid velocity");
		return true;
	}
}
static inline int nint( float val ) { return (int)floor(val + 0.5); }
void EQKAnalyzer::LoadData() {
	// check input params
	// do waveform fitting if model and saclist are input
	_usewaveform = ( !fsaclistR.empty() && !fmodelR.empty() ) || ( !fsaclistL.empty() && !fmodelL.empty() );

	if( _usewaveform ) {	// read in sacs
		// period band
		if( f2<0. || f3<0. || f2>=f3 )
			throw ErrorEA::BadParam(FuncName, "invalid freq range: "+std::to_string(f2)+" - "+std::to_string(f3));

		// lambda function that loads a single sac file into both sac3V and synG
		auto LoadSac = [&]( const std::string& sacname, SynGenerator& synG, std::vector<SacRec3>& sac3V ) {
			SacRec sac(sacname); sac.Load();
			auto& shd = sac.shd;
			if( shd.depmax!=shd.depmax ) return;
			// check/define event location in sac header
			if( shd.evlo<-180 || shd.evla<-90 ) {
				shd.evlo = initlon;
				shd.evla = initlat;
			}
			// check distance
			float dis = sac.Dis();
			if( dis<DISMIN || dis>DISMAX ) return;
			// SNR
			shd.user1 = sac.SNR(dis*0.2, dis*0.5, dis*0.5+500., dis*0.5+1000.);
			if( shd.user1 < SNRMIN ) return;
			// filter
			sac.Resample();	// sample grid alignment
			if( sacRtype == 1 ) sac.Integrate();
			//sac.Filter(f1,f2,f3,f4);
			sac.BandpassCOSFilt(f1,f2,f3,f4);
			// zoom in to the surface wave window
			float tb, te;
			if( dis < 300. ) {
				tb = std::max(dis*0.35-45., 0.);
				te = tb + 90.;
			} else {
				tb = dis*0.2; te = dis*0.5;
			}
			sac.cut(tb, te);
			shd.user2 = tb; shd.user3 = te;
			// FFT
			SacRec sac_am, sac_ph;
			sac.ToAmPh(sac_am, sac_ph);
			// take the amplitude from IFFT for envelope
			SacRec sacEnv; sacEnv.shd = shd;
			sacEnv.FromAmPh(sac_am, sac_ph, 2); 
			shd.user4 = sacEnv.Tpeak();
			// save sac files and put a station record into synGR
			synG.PushbackSta( sac );
			sac3V.push_back( SacRec3{ std::move(sac), std::move(sac_am), std::move(sac_ph) } );
		};

		// Rayleigh
		_synGR.Initialize(fmodelR, fRphvname, fReigname, 'R', 0);
		std::ifstream fin( fsaclistR );
		if( ! fin )
	      throw ErrorEA::BadFile(FuncName, fsaclistR);
		_sac3VR.clear(); _synGR.ClearSta();
		for( std::string line; std::getline(fin, line); ) {
			std::stringstream ss(line); ss >> line;
			LoadSac( line, _synGR, _sac3VR );
		}
		fin.close(); fin.clear();
		float pseudo_per = nint(1./f3) + 0.001*nint(1./f2);
		_dataR.push_back( SDContainer{pseudo_per, R, false} );	// waveform data container
		// Uncertainties
		auto Ispi = spiRM.find(-1.);
		if( Ispi != spiRM.end() ) {
			auto &spi = Ispi->second; auto &sigmaS = _dataR.back().sigmaS;
			sigmaS = AziData{0., spi.sigmaG, spi.sigmaP, spi.sigmaA}; sigmaS = sigmaS * sigmaS;
		}

		// Love
		_synGL.Initialize(fmodelL, fLphvname, fLeigname, 'L', 0);
		fin.open( fsaclistL );
		if( ! fin )
			throw ErrorEA::BadFile(FuncName, fsaclistL);
		_sac3VL.clear(); _synGL.ClearSta();
		for( std::string line; std::getline(fin, line); ) {
			std::stringstream ss(line); ss >> line;
			LoadSac( line, _synGL, _sac3VL );
		}
		_dataL.push_back( SDContainer{pseudo_per, L, false} );	// waveform data container
		// Uncertainties
		Ispi = spiLM.find(-1.);
		if( Ispi != spiLM.end() ) {
			auto &spi = Ispi->second; auto &sigmaS = _dataL.back().sigmaS;
			sigmaS = AziData{0., spi.sigmaG, spi.sigmaP, spi.sigmaA}; sigmaS = sigmaS * sigmaS;
		}

		std::cout<<"### "<<_sac3VR.size()<<"(Rayl) + "<<_sac3VL.size()<<"(Love) sac file(s) loaded. ###"<<std::endl;
	} else {	// read DISP measurements
		auto LoadSDData = [&]( const std::map<float, SinglePeriodInfo>& spiM, std::vector<SDContainer>& dataV, const Dtype t ) {
			dataV.clear(); dataV.reserve( spiM.size() );
			for( const auto& psPair : spiM ) {
				const float per = psPair.first;
				const auto& spi = psPair.second;
				if( spi.fmeasure.empty() || spi.fmapG.empty() || spi.fmapP.empty() ) {
					WarningEA::BadParam(FuncName, "invalid SinglePeriodInfo, skipped.");
					continue;
				}
				// farray[0]: measurement file
				// farray[1] & [2]: either G&P vel_map files or G&P velocities (float for 1D model)
				// farray[3]: list of stations to be used
				float velG, velP;
				SDContainer sd;
				if( FilenameToVel(spi.fmapG, velG) && FilenameToVel(spi.fmapP, velP) ) {
					sd = SDContainer(per, t, true, spi.fmeasure, velG, velP, spi.fstalst, initlon, initlat, DISMIN, DISMAX);	// FTAN measurement container
				} else {
					sd = SDContainer(per, t, true, spi.fmeasure, spi.fmapG, spi.fmapP, spi.fstalst, initlon, initlat, DISMIN, DISMAX);
				}
				if( spi.sigmaG>0 && spi.sigmaP>0 && spi.sigmaA>0 ) {
					sd.sigmaS = AziData{0., spi.sigmaG, spi.sigmaP, spi.sigmaA}; sd.sigmaS = sd.sigmaS * sd.sigmaS;
				}
				dataV.push_back( std::move(sd) );
			}
		};
		// Rayleigh
		LoadSDData( spiRM, _dataR, R );
		// Love
		LoadSDData( spiLM, _dataL, L );

		std::cout<<"### "<<_dataR.size()<<"(Rayl) + "<<_dataL.size()<<"(Love) data periods (measurements+maps) loaded. ###"<<std::endl;
	}

}


/* -------------------- fill the rpR and rpL objects with the current model state -------------------- */
inline std::vector<float> EQKAnalyzer::perRlst() const { return perlst(R); }
inline std::vector<float> EQKAnalyzer::perLlst() const { return perlst(L); }
inline std::vector<float> EQKAnalyzer::perlst(const Dtype dtype) const {
	std::vector<float> perV;
	if( _usewaveform ) {
		perV.push_back(1./f3);
		perV.push_back(1./f2);
	} else {
		const auto &dataV = dtype==R ? _dataR : _dataL;
		for( const auto &data : dataV ) perV.push_back( data.per );
	}
	return perV;
}

// initialize the Analyzer by pre- predicting radpatterns and updating pathpred for all SDContainer based on the current model info
// shift by T multiples according to lower and upper bound. Results not guranteed to be in the range
inline float EQKAnalyzer::ShiftInto( float val, float lb, float ub, float T) const {
	while(val >= ub) val -= T;
	while(val < lb) val += T;
	return val;
}
inline float EQKAnalyzer::BoundInto( float val, float lb, float ub ) const {
	if( val < lb ) val = lb;
	else if( val > ub ) val = ub;
	return val;
}
/*
void EQKAnalyzer::UpdatePredsM( const ModelInfo& minfo, bool updateSource ) {
	UpdatePredsM( minfo, _rpR, _rpL, _dataR, _dataL, _AfactorR, _AfactorL, _source_updated, updateSource );
}
*/
void EQKAnalyzer::UpdatePredsM( const ModelInfo& minfo,	std::vector<SDContainer>& dataR, 
										  std::vector<SDContainer>& dataL, bool updateSource ) const {
	auto rpR = _rpR, rpL = _rpL;
	float AfactorR = _AfactorR, AfactorL = _AfactorL;
	bool source_updated = _source_updated;
	UpdatePredsM( minfo, rpR, rpL, dataR, dataL, AfactorR, AfactorL, source_updated, updateSource );
}
void EQKAnalyzer::UpdatePredsM( const ModelInfo& minfo, RadPattern& rpR, RadPattern& rpL,
										  std::vector<SDContainer>& dataR, std::vector<SDContainer>& dataL,
										  float& AfactorR, float& AfactorL, bool& source_updated, bool updateSource ) const {
	// radpattern
	bool model_updated = false;
	float stk = minfo.stk, dip = minfo.dip, rak = minfo.rak, dep = minfo.dep, M0 = minfo.M0;
	stk = ShiftInto( stk, 0., 360., 360. );	
	dip = BoundInto( dip, 0., 90. );
	rak = ShiftInto( rak, -180., 180., 360. ); dep = BoundInto( dep, 0., 60. );
	model_updated |= rpR.Predict( 'R', fReigname, fRphvname, stk, dip, rak, dep, M0, perRlst() );
	model_updated |= rpL.Predict( 'L', fLeigname, fLphvname, stk, dip, rak, dep, M0, perLlst() );

	if( _usewaveform ) return;

	// SDContainer Tpaths
	for( auto& sdc : dataR )
		model_updated |= sdc.UpdatePathPred( minfo.lon, minfo.lat, minfo.t0 );
	for( auto& sdc : dataL )
		model_updated |= sdc.UpdatePathPred( minfo.lon, minfo.lat, minfo.t0 );

	if( model_updated ) source_updated = false;	// source terms need to be updated
	else if( source_updated ) return;				//	model state didn't change since the last source term update -- return
	if( ! updateSource ) return;						// do not update source for now
	source_updated = true;								// do update source (now)

	// lambda function, works on a single SDContainer vector
	auto Usource = [&]( std::vector<SDContainer>& dataV, float& Afactor ) {
		// update source predictions
		for( auto& sdc : dataV ) {
			const auto& rp = sdc.type==R ? rpR : rpL;
			sdc.UpdateSourcePred( rp );
		}
		/*	Amplitudes are now correctly scaled by RadPattern::GetPred
		// rescale source amplitudes to match the data
		std::vector<float> ampratioV;	
		ampratioV.reserve( dataV.size() * dataV.at(0).size() );
		for( auto& sdc : dataV ) sdc.ComputeAmpRatios( ampratioV );
		Afactor = std::accumulate(ampratioV.begin(), ampratioV.end(), 0.) / ampratioV.size();
		//std::cerr<<"   Afactor for Rayleigh: "<<Afactor<<" (was 50000.)\n";
		for( auto& sdc : dataV )	// scale source amplitudes to match the observations
			sdc.AmplifySource( Afactor );
		*/
	};

	// SDContainer R: source terms
	if( dataR.size() > 0 ) Usource( dataR, AfactorR );

	// SDContainers L: source terms
	if( dataL.size() > 0 ) Usource( dataL, AfactorL );
}


/* -------------------- compute the total chi-square misfit based on the current data state and the input model info -------------------- */
void EQKAnalyzer::chiSquareM( ModelInfo minfo, float& chiS, int& N ) const {
	// check/correct input params 
	minfo.Correct();

	int Rsize = _dataR.size(), Lsize = _dataL.size();
	if( Rsize==0 && Lsize==0 )
		throw ErrorEA::EmptyData(FuncName, "Rsize && Lsize");
	if( datatype == Undefined )
		throw ErrorEA::BadParam(FuncName, "Undefined datatype");

	// data flags
	bool useG = _useG, useP = _useP, useA = _useA;
	bool RFlag = (datatype==B || datatype==R);
	bool LFlag = (datatype==B || datatype==L);
	if( _isInit ) {
		useG = true; useA = useP = false;
		if( LFlag ) { RFlag = false; }
	}

	// initialize ADAdder
	ADAdder adder( useG, useP, useA );
	const int Nadd = useG + useP + useA;

	// make a copy of the data and RadPattern objs (for multi-threading),
	// update both path and source predictions for all SDcontainers
	auto dataR = _dataR, dataL = _dataL;
	UpdatePredsM( minfo, dataR, dataL, true );

	// lambda: computes chi-square of a single wavetype (with a vector of SDContainer)
	chiS = 0.; N = 0; //wSum = 0.;
	auto chiSM = [&](std::vector<SDContainer> &SDV) {
		for( auto& sdc : SDV ) {
			// compute misfits on all stations,
			// average in each (20 degree) bin,
			// and store the results into AziData vectors
			std::vector<AziData> adVmean, adVvar;
			auto getVar = [&](const int i){ return dynamicVar ? adVvar[i] : sdc.sigmaS; };
			//std::cerr<<"debug(Check Sigmas): "<<dynamicVar<<" "<<getVar(0)<<" "<<getVar(10)<<std::endl;
			sdc.BinAverage( adVmean, adVvar, true, true, dynamicVar );	// correct for 2pi, data taken as FTAN measurements
			// add to chi-square misfit
			for( int i=0; i<adVmean.size(); i++ ) {
				const auto& admean = adVmean[i];
				//const auto& advar = adVvar[i];
				// (mis * mis / variance)
				chiS += adder( (admean * admean)/getVar(i) );
				N += Nadd;
			}
			//std::cerr<<"Rayleigh at per="<<sdc.per<<": "<<chiS<<" "<<N<<" "<<chiS/(N-8.)<<"\n";
		}
	};

	// accumulate for Rayleigh
	if( Rsize>0 && RFlag ) chiSM(dataR);

	// accumulate for Love
	if( Lsize>0 && LFlag ) chiSM(dataL);

}


// given an observed waveform (sac, sac_am, sac_ph), generate synthetic and compute misfit
SacRec EQKAnalyzer::ComputeSyn(const SacRec &sac, SynGenerator &synG) const {
	const auto &shd = sac.shd;
	SacRec sacSZ, sacSR, sacST;
	int nptsS = ceil( (shd.user3-synG.minfo.t0)/shd.delta ) + 1;
	bool synsuc = synG.ComputeSyn( sac.stname(), shd.stlo, shd.stla, nptsS, shd.delta, 
											 sacSZ, sacSR, sacST, rotateSyn, f1, f2, f3, f4 );
	Dtype type = synG.type=='R' ? R : L;
	if( ! synsuc ) {
		(type==R ? sacSZ : sacST).Write("debug_Syn_"+sac.stname()+".SAC");
		throw ErrorEA::InternalException(FuncName, "bad sac from ComputeSyn for station "+sac.stname());
		//WarningEA::Other(FuncName, "failed to produce/invalid synthetics for station "+sac.stname() );
		//return;
	}

	// clear-up/select sac
	/*
	if( type == R ) {
		sacSR.clear(); sacST.clear();	// clear unnecessory channels
	} else {
		sacSZ.clear(); sacSR.clear();
	}
	*/

	// check for bad sac (done already in ComputeSyn!)
	//if( sacS.shd.depmax!=sacS.shd.depmax ) return;

	return type==R ? sacSZ : sacST;
}

StaData EQKAnalyzer::WaveformMisfit( const SacRec3 &sac3, SynGenerator &synG ) const {
	// references to data sacs
	const SacRec &sacM = sac3[0], &sac_am1 = sac3[1], &sac_ph1 = sac3[2];

	// produce synthetic
	SacRec sacS = ComputeSyn(sacM, synG);

	sacS.Resample();	// important! shift to regular sampling grids

	// zoom in to the surface wave window
	auto& shdM = sacM.shd;
	//float tb = shdM.user2, te = shdM.user3;
	sacS.cut( shdM.user2, shdM.user3 );

	// FFT into freq-domain
	SacRec sac_am2, sac_ph2;
	sacS.ToAmPh(sac_am2, sac_ph2);
	// take the amplitude from IFFT for envelope
	sacS.FromAmPh(sac_am2, sac_ph2, 2); 

	// group time shift
	float grTShift = shdM.user4 - sacS.Tpeak();

	// time-domain correlation
	//std::cout<<shdM.stlo<<" "<<shdM.stla<<" "<<sacM.Correlation( sacS, tmin, tmax )<<"   "<<sacM.stname()<<" "<<sacS.stname()<<" "<<sacS.Dis()<<" "<<sacS.Azi()<<" CC_sigtime "<<sacM.fname<<"\n";

	// freq-domain rms
	float AmpTheory; sac_am2.Mean(f2, f3, AmpTheory);
	sac_am2.Subf(sac_am1); sac_ph2.Subf(sac_ph1);
	// correct for 2pi
	const float TWO_PI = M_PI * 2.;
	auto correct2PI = [&](float& val) {
		if( val >= M_PI ) val -= TWO_PI;
		else if( val < -M_PI ) val += TWO_PI;
	};
	sac_ph2.Transform( correct2PI );
	float rms_am = -sac_am2.RMSAvg(f2, f3, true), rms_ph = -sac_ph2.RMSAvg(f2, f3, true);	// signed rms
	//std::cerr<<grTShift<<" "<<rms_ph<<" "<<rms_am<<"   "<<rms_am+AmpTheory<<" "<<AmpTheory<<std::endl;
	//amp_sum += rms_am; pha_sum += rms_ph; N++;
	// store rms misfits as StaData. Put amp of synthetic in .Asource for computing variance later
	StaData sd(sacS.Azi(), shdM.stlo, shdM.stla, grTShift, rms_ph, rms_am+AmpTheory, 0.); sd.Asource = AmpTheory;
	return sd;
	//data.push_back( sd );
	//std::cout<<"RMS_amp = "<<rms_am<<"   RMS_pha = "<<rms_ph<<std::endl;
	//std::cerr<<lon<<" "<<lat<<" "<<rms_ph<<"   "<<sacM.stname()<<" "<<sacS.stname()<<" "<<sacS.Dis()<<" "<<sacS.Azi()<<" rms_pha "<<sacM.fname<<"\n";
	//std::cerr<<lon<<" "<<lat<<" "<<100.*rms_am/AmpTheory<<"   "<<sacM.stname()<<" "<<sacS.stname()<<" "<<sacS.Dis()<<" "<<sacS.Azi()<<" rms_amp "<<sacM.fname<<"\n";
}

void EQKAnalyzer::UpdatePredsW( const ModelInfo& minfo, std::vector<SDContainer> &dataR, std::vector<SDContainer> &dataL ) const {
	// lambda: update misfits for a single wavetype
	auto updatePreds = [&]( SynGenerator synG, const std::vector<SacRec3> &sac3V, SDContainer &data ) {
		// prepare SynGenerator
		//auto synGR = _synGR;
		synG.SetEvent( minfo );
		//float pseudo_per = nint(1./f3) + 0.001*nint(1./f2);
		//SDContainer data( pseudo_per, synG.type=='R' ? R : L, false );	// waveform data container
		for( const auto& sac3 : sac3V ) data.push_back( WaveformMisfit(sac3, synG) );
		//std::cout<<"average misfits = "<<sqrt(amp_sum/(N-1))<<" "<<sqrt(pha_sum/(N-1))<<std::endl;
		data.Sort(); data.UpdateAziDis( minfo.lon, minfo.lat );
	};
	
	if( ! dataR.empty() ) updatePreds( _synGR, _sac3VR, dataR[0] );
	if( ! dataL.empty() ) updatePreds( _synGL, _sac3VL, dataL[0] );
}

void EQKAnalyzer::chiSquareW( ModelInfo minfo, float& chiS, int& N ) const {
	// check/correct input params 
	minfo.Correct();

	int Rsize = _dataR.size(), Lsize = _dataL.size();
	if( Rsize!=1 && Lsize!=1 ) throw ErrorEA::EmptyData(FuncName, "Rsize && Lsize");
	if( datatype == Undefined ) throw ErrorEA::BadParam(FuncName, "Undefined datatype");

	// data flags
	bool RFlag = (datatype==B || datatype==R);
	bool LFlag = (datatype==B || datatype==L);

	// initialize ADAdder
	bool useG = _useG, useP = _useP, useA = _useA;
	if( _isInit ) {
		useG = true; useA = useP = false;
		//if( LFlag ) { RFlag = false; }
	}
	ADAdder adder( useG, useP, useA );
	const int Nadd = useG + useP + useA;

	// make a copy of the data (for multi-threading),
	// update predictions (into misfits) for all SDcontainers
	std::vector<SDContainer> dataR, dataL;
	if(! _dataR.empty()) { dataR.push_back(_dataR[0]); dataR.back().clear(); } // TODO: a SDContainer::method to copy header info only
	if(! _dataL.empty()) { dataL.push_back(_dataL[0]); dataL.back().clear(); }
	UpdatePredsW( minfo, dataR, dataL );

	// lambda: computes chi square for a single wavetype
	auto chiSW = [&]( SynGenerator synG, SDContainer& data ) {
/*
		// prepare SynGenerator
		//auto synGR = _synGR;
		synG.SetEvent( minfo );
		float pseudo_per = nint(1./f3) + 0.001*nint(1./f2);
		SDContainer data( pseudo_per, synG.type=='R' ? R : L, false );	// waveform data container

		//float amp_sum = 0., pha_sum = 0.;
		for( const auto& sac3 : sac3V ) data.push_back( WaveformMisfit(sac3, synG) );
		//std::cout<<"average misfits = "<<sqrt(amp_sum/(N-1))<<" "<<sqrt(pha_sum/(N-1))<<std::endl;
		data.Sort();
*/

		// compute chi square
		std::vector<AziData> adVmean, adVvar;
		auto getVar = [&](const int i){ return dynamicVar ? adVvar[i] : data.sigmaS; };
		data.BinAverage( adVmean, adVvar, false, false, dynamicVar );	// do not correct for 2 pi, data taken as waveform rms misfits
		for( int i=0; i<adVmean.size(); i++ ) {
			const auto& admean = adVmean[i];
			//const auto& advar = adVvar[i];
			//std::cerr<<"admean = "<<admean<<"   advar = "<<advar<<"   adder = "<<(admean*admean)/advar<<std::endl;
			// (mis * mis / variance)
			chiS += adder( (admean * admean)/getVar(i) );
			N += Nadd;
		}
/*
		// feed dataout/dataout if requested
		if( filldata ) {
			data.UpdateAziDis(minfo.lon, minfo.lat);
			dataout = std::move(data);
		}
*/
	};

	// compute Rayleigh and Love
	chiS = 0.; N = 0;
	if( Rsize==1 && RFlag ) chiSW( _synGR, dataR[0] );
	if( Lsize==1 && LFlag ) chiSW( _synGL, dataL[0] );

}


// output real (processed) and synthetic waveforms when the waveform fitting method is used
void EQKAnalyzer::OutputWaveforms( const ModelInfo& minfo, const std::string& outdir ) {
	if( ! _usewaveform ) return;
		//throw ErrorEA::BadParam(FuncName, "non waveform-fitting");
	if( outdir.empty() )
		throw ErrorEA::BadParam(FuncName, "empty outsac_dir");
	// data flags
	bool RFlag = (datatype==B || datatype==R);
	bool LFlag = (datatype==B || datatype==L);

   MKDirs( outdir );

	// lambda: output all waveforms for a single wavetype
	auto OutputW = [&]( SynGenerator& synG, std::vector<SacRec3>& sac3V ) {
		// prepare SynGenerator
		//auto& synGR = _synGR;
		Dtype type = synG.type=='R' ? R : L;
		synG.SetEvent( minfo );

		for( auto& sac3 : sac3V ) {
			auto &sacM = sac3[0]; auto &shdM = sacM.shd;
			// produce synthetic
			SacRec sacS = ComputeSyn(sacM, synG);
			sacS.Resample();	// shift to regular sampling grids
			sacS.cut( shdM.user2, shdM.user3 );

			std::string sacname = outdir + "/" + sacS.stname() + "." + sacS.chname();
			sacM.Write( sacname + "_real.SAC" ); sacS.Write( sacname + "_syn.SAC" );
		}
	};

	// call lambda for Rayleigh and Love
	if( RFlag ) OutputW( _synGR, _sac3VR );
	if( LFlag ) OutputW( _synGL, _sac3VL );

}


// output variances for R & L all periods into a single file
void EQKAnalyzer::EstimateSigmas( const ModelInfo& minfo ) {
	// prepare/update data for the input minfo
	if( _usewaveform ) {	// fill data vectors for the waveform fitting method
		UpdatePredsW(minfo);
	} else {	// update both path and source predictions for all SDcontainers
		UpdatePredsM( minfo, true );
	}

	for( auto& sdc : _dataR ) sdc.ComputeVar();

	for( auto& sdc : _dataL ) sdc.ComputeVar();
}
void EQKAnalyzer::OutputSigmas() const {
	// open output file
	std::ofstream fout( outname_sigmas, std::ofstream::app );
	//fout<<"# [ minfo = "<<minfo<<" ]\n";
	fout<<"Rayl:\n";
	for( auto& sdc : _dataR ) fout<<sqrt(sdc.sigmaS)<<"\n";
	fout<<"Love:\n";
	for( auto& sdc : _dataL ) fout<<sqrt(sdc.sigmaS)<<"\n";
	fout<<"\n";
}


/* -------------------- Output the data and predictions based on the input model info -------------------- */
void EQKAnalyzer::OutputFits( ModelInfo minfo ) {
	// check/correct input params
	minfo.Correct();

	// prepare/update data for the input minfo
	if( _usewaveform ) {	// fill data vectors for the waveform fitting method
		UpdatePredsW(minfo);
	} else {	// update both path and source predictions for all SDcontainers
		UpdatePredsM( minfo, true );
	}
	bool isFTAN = ! _usewaveform;

	// check data sizes
	int Rsize = _dataR.size(), Lsize = _dataL.size();
	if( Rsize==0 && Lsize==0 )
		throw ErrorEA::EmptyData(FuncName, "Rsize && Lsize");

	// lambda function for outputing fit
	auto outputF = [&]( SDContainer& sdc ) {
		const float per = sdc.per;
		// choose outlist by Dtype
		const auto& outlist = sdc.type==R ? outlist_RF : outlist_LF;
		const auto& rp = sdc.type==R ? _rpR : _rpL;
		const float Afactor = sdc.type==R ? _AfactorR : _AfactorL;
		std::string outname;
		if( _usewaveform ) {
			if( outlist.size() == 0 ) return;
			outname = outlist.begin()->second;
		} else { // outname for the current period
			if( outlist.find(per) == outlist.end() ) return;
			outname = outlist.at(per);
		}
		// average in each (20 degree) bin,
		std::vector<AziData> adVmean, adVvar;
		auto getVar = [&](const int i){ return dynamicVar ? adVvar[i] : sdc.sigmaS; };
		sdc.BinAverage( adVmean, adVvar, isFTAN, isFTAN, dynamicVar );	// 2pi gets corrected here
		// output data and predictions at each station (append at the end)
		std::ofstream foutsta( outname + "_sta", std::ofstream::app );
		foutsta<<"# [ minfo = "<<minfo<<" ]\n";
		sdc.PrintAll( foutsta, true );	// normalize amplitude while printing
		foutsta << "\n\n";
		// output misfits and source preds at each bin azi (append at the end)
		std::ofstream foutbin( outname + "_bin", std::ofstream::app );
		foutbin<<"# [ minfo = "<<minfo<<" ]\n";
		for( int i=0; i<adVmean.size(); i++ ) {
			const auto& admean = adVmean[i];
			//const auto adstd = sqrt( adVvar[i] );
			const auto adstd = sqrt( getVar(i) );
			float Gsource, Psource, Asource;
			if( _usewaveform ) {
				Gsource = Psource = Asource = 0.;
			} else {
				rp.GetPred( per, admean.azi,	Gsource, Psource, Asource );
							//	1000., 0., rp.cAmp(per)[0], admean.Gdata );
				//Asource *= Afactor;
			}
			// averaged-normalized amplitudes are stored in admean.user!
			foutbin<<admean.azi<<"  "<<Gsource<<" "<<admean.Gdata<<" "<<adstd.Gdata
				<<"  "<<Psource<<" "<<admean.Pdata<<" "<<adstd.Pdata
				<<"  "<<admean.user<<" "<<admean.Adata<<" "<<adstd.Adata<<"\n";
		}
		foutbin << "\n\n";
	};

	// main loop for Rayleigh
	if( Rsize>0 ) for( auto& sdc : _dataR ) outputF(sdc);

	// main loop for Love
	if( Lsize>0 ) for( auto& sdc : _dataL ) outputF(sdc);

}


/* -------------------- compute and output misfits, separately, for group, phase, and amplitudes -------------------- */
void EQKAnalyzer::OutputMisfits( ModelInfo minfo ) {
	// check/correct input params
	minfo.Correct();

	// prepare/update data for the input minfo
	if( _usewaveform ) {	// fill data vectors for the waveform fitting method
		UpdatePredsW(minfo);
	} else {	// update both path and source predictions for all SDcontainers
		UpdatePredsM( minfo, true );
	}
	bool isFTAN = ! _usewaveform;

	// check data sizes
	int Rsize = _dataR.size(), Lsize = _dataL.size();
	if( Rsize==0 && Lsize==0 )
		throw ErrorEA::EmptyData(FuncName, "Rsize && Lsize");

	// prepare output file
	std::ofstream fout( outname_misAll, std::ofstream::app );
	fout<<"# [ minfo = "<<minfo<<" ]\n";
	fout<<"wtype per  rmsG L1G rchisG  rmsP L1P rchisP  rmsA L1A rchisA\n";

	// lambda function for outputing misfits
	auto outputM = [&]( SDContainer& sdc ) {
		const float per = sdc.per;
		// average in each (20 degree) bin,
		std::vector<AziData> adVmean, adVvar;
		auto getVar = [&](const int i){ return dynamicVar ? adVvar[i] : sdc.sigmaS; };
		sdc.BinAverage( adVmean, adVvar, isFTAN, isFTAN );
		// output misfits and source preds at each bin azi (append at the end)
		AziData misL1{0.}, misL2{0.}, chiS{0.}; 
		float nbin = adVmean.size();
		for( int i=0; i<nbin; i++ ) {
			const auto& admean = adVmean[i];
			misL1 += fabs(admean); misL2 += admean * admean;
			chiS += (admean*admean) / getVar(i);
		}
		misL1 /= nbin; misL2 = sqrt( misL2/(nbin-1.) );
		chiS /= nbin;
		char type = sdc.type==R ? 'R' : 'L';
		fout << type << " " << per << "  "
			<< misL2.Gdata << " " << misL1.Gdata << " " <<  chiS.Gdata << "  "
			<< misL2.Pdata << " " << misL1.Pdata << " " <<  chiS.Pdata << "  "
			<< misL2.Adata << " " << misL1.Adata << " " <<  chiS.Adata << "\n";
	};

	// main loop for Rayleigh
	if( Rsize>0 ) for( auto& sdc : _dataR ) outputM(sdc);

	// main loop for Love
	if( Lsize>0 ) for( auto& sdc : _dataL ) outputM(sdc);

}


/* -------------------- output source predictions (continuously in azimuth, for group, phase, and amplitudes) into single file for R/L waves -------------------- */
void EQKAnalyzer::OutputSourcePatterns( const ModelInfo& mi ) {
	if( outname_srcR.empty() && outname_srcL.empty() ) return;
	// call UpdatePredsM for Afactors (needs a new, less redundant method here!)
	UpdatePredsM( mi, true );
	/*
	float stk = mi.stk, dip = mi.dip, rak = mi.rak, dep = mi.dep;
	//stk = ShiftInto( stk, 0., 360., 360. );	dip = BoundInto( dip, 0., 90. );
	//rak = ShiftInto( rak, -180., 180., 360. ); dep = BoundInto( dep, 0., 60. );
	_rpR.Predict( 'R', fReigname, fRphvname, stk, dip, rak, dep, perRlst() );
	_rpL.Predict( 'L', fLeigname, fLphvname, stk, dip, rak, dep, perLlst() );
	*/

	// output source predictions (R)
	if( ! outname_srcR.empty() ) _rpR.OutputPreds(outname_srcR, _AfactorR);
	if( ! outname_srcL.empty() ) _rpL.OutputPreds(outname_srcL, _AfactorL);
}

