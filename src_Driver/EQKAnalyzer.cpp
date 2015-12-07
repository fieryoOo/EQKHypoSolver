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
	if( stmp == "fRse" ) succeed = buff >> fReigname;
	else if( stmp == "fRsp" ) succeed = buff >> fRphvname;
	else if( stmp == "fLse" ) succeed = buff >> fLeigname;
	else if( stmp == "fLsp" ) succeed = buff >> fLphvname;
	else if( stmp == "weightR_Loc" ) succeed = buff >> weightR_Loc;
	else if( stmp == "weightL_Loc" ) succeed = buff >> weightL_Loc;
	else if( stmp == "weightR_Foc" ) succeed = buff >> weightR_Foc;
	else if( stmp == "weightL_Foc" ) succeed = buff >> weightL_Foc;
	else if( stmp == "noG" ) { succeed = true; _useG = false; }
	else if( stmp == "noP" ) { succeed = true; _useP = false; }
	else if( stmp == "noA" ) { succeed = true; _useA = false; }
	else if( stmp == "SNRMIN" ) succeed = buff >> SNRMIN;
	else if( stmp == "DISMIN" ) succeed = buff >> DISMIN;
	else if( stmp == "DISMAX" ) succeed = buff >> DISMAX;
	else if( stmp == "lon" ) { 
		succeed = buff >> initlon; 
		if( succeed && initlon<0.) initlon += 360.; 
	}
	else if( stmp == "lat" ) succeed = buff >> initlat;
	else if( stmp == "fmodelR" ) succeed = buff >> fmodelR;
	else if( stmp == "fmodelL" ) succeed = buff >> fmodelL;
	else if( stmp == "fsaclistR" ) {
		succeed = buff >> fsaclistR >> sacRtype;
		if( sacRtype!=0 && sacRtype!=1 )
			throw ErrorEA::BadParam( FuncName, "sacRtype: expecting either 0 (displacement) or 1 (velocity)");
	}
	else if( stmp == "fsaclistL" ) {
		succeed = buff >> fsaclistL >> sacLtype;
		if( sacLtype!=0 && sacLtype!=1 )
			throw ErrorEA::BadParam( FuncName, "sacLtype: expecting either 0 (displacement) or 1 (velocity)");
	}
	else if( stmp == "permin" ) {
		float permin;
		succeed = buff >> permin;
		f3 = 1./permin; f4 = 1.1/permin;
	}
	else if( stmp == "permax" ) {
		float permax;
		succeed = buff >> permax;
		f1 = 0.9/permax; f2 = 1./permax;
	}
	else if( stmp == "indep" ) { 
		succeed = buff >> _indep_factor;
		if( _indep_factor<=0. || _indep_factor>1. )
			throw ErrorEA::BadParam( FuncName, "indep factor out of range: "+std::to_string(_indep_factor) );
	}
	else if( stmp == "dflag" ) {
		succeed = buff >> datatype_name;
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
		succeed = buff >> stmp1 >> stmp2 >> stmp3 >> per;	// fmeasure, fmapG, fmapP, per
		buff >> stmp4;													// station list (optional)
		if( succeed ) {
			if( fRlist.find(per) != fRlist.end() ) return -3;	// period already exists
			fRlist[per] = std::array<FileName, 4>{stmp1, stmp2, stmp3, stmp4};
		}
	}
	else if( stmp == "fLm" ) {
		FileName stmp1, stmp2, stmp3, stmp4;
		float per;
		succeed = buff >> stmp1 >> stmp2 >> stmp3 >> per;	// fmeasure, fmapG, fmapP, per
		buff >> stmp4;													// station list (optional)
		if( succeed ) {
			if( fLlist.find(per) != fLlist.end() ) return -3;	// period already exists
			fLlist[per] = std::array<FileName, 4>{stmp1, stmp2, stmp3, stmp4};
		}
	}
	else if( stmp == "fmisL" ) {
		FileName& outname = outname_misL;
		succeed = buff >> outname;
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
		succeed = buff >> outname;
		if( succeed && MoveF ) outname.SaveOld();
	}
	else if( stmp == "fmisAll" ) {
		FileName& outname = outname_misAll;
		succeed = buff >> outname;
		if( succeed && MoveF ) outname.SaveOld();
	}
	else if( stmp == "fpos" ) {
		FileName& outname = outname_pos;
		succeed = buff >> outname;
		if( succeed && MoveF ) outname.SaveOld();
	}
	else if( stmp == "fsrcR" ) {
		FileName& outname = outname_srcR;
		succeed = buff >> outname;
		if( succeed && MoveF ) outname.SaveOld();
	}
	else if( stmp == "fsrcL" ) {
		FileName& outname = outname_srcL;
		succeed = buff >> outname;
		if( succeed && MoveF ) outname.SaveOld();
	}
	else if( stmp == "dirsac" ) {
		FileName& dirname = outdir_sac;
		succeed = buff >> dirname;
		if( succeed && MoveF ) dirname.SaveOld();
	}
	else if( stmp == "ffitR" ) {
		FileName outname;
		float per;
		succeed = buff >> outname >> per;
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
		succeed = buff >> outname >> per;
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
/*
	else if( stmp == "fpredR" ) {
		FileName outname;
		float per;
		succeed = buff >> outname >> per;
		if( succeed ) {
			outlist_RP[per] = outname;
			//if( MoveF ) outname.SaveOld();
		}
	}
	else if( stmp == "fpredL" ) {
		FileName outname;
		float per;
		succeed = buff >> outname >> per;
		if( succeed ) {
			outlist_LP[per] = outname;
			//if( MoveF ) outname.SaveOld();
		}
	} */
	else return -1;
	if( succeed ) return 1;
	return 0;
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
void EQKAnalyzer::LoadData() {
	// check input params
	// do waveform fitting if model and saclist are input
	_usewaveform = ( !fsaclistR.empty() && !fmodelR.empty() ) || ( !fsaclistL.empty() && !fmodelL.empty() );

	if( _usewaveform ) {	// read in sacs
		// period band
		if( f2<0. || f3<0. || f2>=f3 )
			throw ErrorSC::BadParam(FuncName, "invalid freq range: "+std::to_string(f2)+" - "+std::to_string(f3));

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
	      throw ErrorSC::BadFile(FuncName, fsaclistR);
		_sac3VR.clear(); _synGR.ClearSta();
		for( std::string line; std::getline(fin, line); ) {
			std::stringstream ss(line); ss >> line;
			LoadSac( line, _synGR, _sac3VR );
		}
		fin.close(); fin.clear();

		// Love
		_synGL.Initialize(fmodelL, fLphvname, fLeigname, 'L', 0);
		fin.open( fsaclistL );
		if( ! fin )
			throw ErrorSC::BadFile(FuncName, fsaclistL);
		_sac3VL.clear(); _synGL.ClearSta();
		for( std::string line; std::getline(fin, line); ) {
			std::stringstream ss(line); ss >> line;
			LoadSac( line, _synGL, _sac3VL );
		}

		std::cout<<"### "<<_sac3VR.size()<<"(Rayl) + "<<_sac3VL.size()<<"(Love) sac file(s) loaded. ###"<<std::endl;
	} else {	// read DISP measurements
		auto LoadSDData = [&]( const std::map<float, std::array<FileName, 4> >& flist, std::vector<SDContainer>& dataV, const Dtype t ) {
			dataV.clear(); dataV.reserve( flist.size() );
			for( const auto& f : flist ) {
				const float per = f.first;
				const auto& farray = f.second;
				// farray[0]: measurement file
				// farray[1] & [2]: either G&P vel_map files or G&P velocities (float for 1D model)
				// farray[3]: list of stations to be used
				float velG, velP;
				SDContainer sd;
				if( FilenameToVel(farray[1], velG) && FilenameToVel(farray[2], velP) ) {
					sd = SDContainer(per, t, farray[0], velG, velP, farray[3], initlon, initlat, DISMIN, DISMAX);
				} else {
					sd = SDContainer(per, t, farray[0], farray[1], farray[2], farray[3], initlon, initlat, DISMIN, DISMAX);
				}
				dataV.push_back( std::move(sd) );
			}
		};
		// Rayleigh
		LoadSDData( fRlist, _dataR, R );
		// Love
		LoadSDData( fLlist, _dataL, L );

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
		auto dataV = dtype==R ? _dataR : _dataL;
		for( const auto& data : dataV ) perV.push_back( data.per );
	}
	return perV;
}

// initialize the Analyzer by pre- predicting radpatterns and updating pathpred for all SDContainer based on the current model info
void EQKAnalyzer::PredictAll( const ModelInfo& minfo, bool updateSource ) {
	PredictAll( minfo, _rpR, _rpL, _dataR, _dataL, _AfactorR, _AfactorL, _source_updated, updateSource );
}
void EQKAnalyzer::PredictAll( const ModelInfo& minfo,	std::vector<SDContainer>& dataR, 
										std::vector<SDContainer>& dataL, bool updateSource ) const {
	auto rpR = _rpR, rpL = _rpL;
	float AfactorR = _AfactorR, AfactorL = _AfactorL;
	bool source_updated = _source_updated;
	PredictAll( minfo, rpR, rpL, dataR, dataL, AfactorR, AfactorL, source_updated, updateSource );
}
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
void EQKAnalyzer::PredictAll( const ModelInfo& minfo, RadPattern& rpR, RadPattern& rpL,
										std::vector<SDContainer>& dataR, std::vector<SDContainer>& dataL,
										float& AfactorR, float& AfactorL, bool& source_updated, bool updateSource ) const {
	// radpattern
	bool model_updated = false;
	float stk = minfo.stk, dip = minfo.dip, rak = minfo.rak, dep = minfo.dep, M0 = minfo.M0;
	stk = ShiftInto( stk, 0., 360., 360. );	dip = BoundInto( dip, 0., 90. );
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
		/*	Amplitudes are now correcly scaled by RadPattern::GetPred
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
	PredictAll( minfo, dataR, dataL, true );

	// main loop for Rayleigh
	chiS = 0.; N = 0; //wSum = 0.;
	if( Rsize>0 && RFlag ) {
		for( auto& sdc : dataR ) {
			// compute misfits on all stations,
			// average in each (20 degree) bin,
			// and store the results into AziData vectors
			std::vector<AziData> adVmean, adVvar;
			sdc.BinAverage( adVmean, adVvar );	// correct for 2pi, data taken as FTAN measurements
			// add to chi-square misfit
			for( int i=0; i<adVmean.size(); i++ ) {
				const auto& admean = adVmean[i];
				const auto& advar = adVvar[i];
				//AziData adtmp = 1./(adstd*adstd);
				//wSum += adder( adtmp ); // (adtmp.Gdata + adtmp.Pdata + adtmp.Adata);
				//adtmp *= (admean * admean);
				// (mis * mis / variance)
				chiS += adder( (admean * admean)/advar );
				N += Nadd;
			}
			//std::cerr<<"Rayleigh at per="<<sdc.per<<": "<<chiS<<" "<<N<<" "<<chiS/(N-8.)<<"\n";
		}
	}

	// main loop for Love
	if( Lsize>0 && LFlag ) {
		for( auto& sdc : dataL ) {
			// compute misfits on all stations,
			// average in each (20 degree) bin,
			// and store the results into AziData vectors
			std::vector<AziData> adVmean, adVvar;
			sdc.BinAverage( adVmean, adVvar );	// correct for 2pi, data taken as FTAN measurements
			// compute chi-square misfit
			for( int i=0; i<adVmean.size(); i++ ) {
				const auto& admean = adVmean[i];
				const auto& advar = adVvar[i];
				// (mis * mis / variance)
				chiS += adder( (admean * admean)/advar );
				N += Nadd;
			}
			//std::cerr<<"Love at per="<<sdc.per<<": "<<chiS<<" "<<N<<" "<<chiS/(N-8.)<<"\n";
		}
	}

}

/* use SacRec::Tpeak instead
float EQKAnalyzer::Tpeak( const SacRec& sac ) const {
	//float tmin, tmax, min, max;
	//sac.MinMax(tb, te, tmin, min, tmax, max);
	//return fabs(min)>fabs(max) ? tmin : tmax;	// for signal
	// for envelope
	int imin, imax;
	sac.MinMax(imin, imax);
	if( imax==0 || imax==sac.shd.npts-1 ) {
		return sac.Time(imax);
	} else {
		float* sigsac = sac.sig.get();
		PointC p1(sac.Time(imax-1), sigsac[imax-1]);
		PointC p2(sac.Time(imax),   sigsac[imax]  );
		PointC p3(sac.Time(imax+1), sigsac[imax+1]);
		return Parabola( p1, p2, p3 ).Vertex().x;
	}
}
*/

static inline int nint( float val ) { return (int)floor(val + 0.5); }
void EQKAnalyzer::chiSquareW( ModelInfo minfo, float& chiS, int& N, bool filldata, SDContainer& dataRout, SDContainer& dataLout ) const {
	// check/correct input params
	minfo.Correct();

	// data flags
	bool RFlag = (datatype==B || datatype==R);
	bool LFlag = (datatype==B || datatype==L);

	// initialize ADAdder
	bool useG = _useG, useP = _useP, useA = _useA;
	if( _isInit ) {
		useG = true; useA = useP = false;
		//if( Lsize > 0 ) {	RFlag = false; LFlag = true; }
	}
	ADAdder adder( useG, useP, useA );
	const int Nadd = useG + useP + useA;


	// lambda: computes chi square for a single wavetype
	auto chiSW = [&]( SynGenerator synG, const std::vector<SacRec3>& sac3V, SDContainer& dataout ) {
		// prepare SynGenerator
		//auto synGR = _synGR;
		synG.SetEvent( minfo );
		Dtype type = synG.type=='R' ? R : L;
		float pseudo_per = nint(1./f3) + 0.001*nint(1./f2);
		SDContainer data( pseudo_per, type );

		//float amp_sum = 0., pha_sum = 0.;
		for( const auto& sac3 : sac3V ) {
			// references to data sacs
			const SacRec &sacM = sac3[0], &sac_am1 = sac3[1], &sac_ph1 = sac3[2];
			auto& shdM = sacM.shd;

			// produce synthetic
			SacRec sacSZ, sacSR, sacST;
			float lon=shdM.stlo, lat=shdM.stla;
			int nptsS = ceil( (shdM.user3-minfo.t0)/shdM.delta ) + 1;
			bool synsuc = synG.ComputeSyn( sacM.stname(), lon, lat, nptsS, shdM.delta, 
													 f1,f2,f3,f4, sacSZ, sacSR, sacST, rotateSyn );
			if( ! synsuc ) {
				WarningEA::Other(FuncName, "failed to produce/invalid synthetics for station "+sacM.stname() );
				continue;
			}

			// clear-up/select sac
			if( type == R ) {
				sacSR.clear(); sacST.clear();	// ignore unnecessory channel
			} else {
				sacSZ.clear(); sacSR.clear();
			}
			SacRec& sacS = type==R ? sacSZ : sacST;

			// check for bad sac 
			if( sacS.shd.depmax!=sacS.shd.depmax ) continue;
			sacS.Resample();	// important! shift to regular sampling grids

			// zoom in to the surface wave window
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
			//std::cout<<lon<<" "<<lat<<" "<<sacM.Correlation( sacS, tmin, tmax )<<"   "<<sacM.stname()<<" "<<sacS.stname()<<" "<<sacS.Dis()<<" "<<sacS.Azi()<<" CC_sigtime "<<sacM.fname<<"\n";

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
			float rms_am = sac_am2.RMSAvg(f2, f3), rms_ph = sac_ph2.RMSAvg(f2, f3);
			//amp_sum += rms_am; pha_sum += rms_ph; N++;
			// store rms misfits as StaData. Put amp of synthetic in .Asource for computing variance later
			StaData sd(sacS.Azi(), lon, lat, grTShift, rms_ph, rms_am+AmpTheory, 0.); sd.Asource = AmpTheory;
			data.push_back( sd );
			//std::cout<<"RMS_amp = "<<rms_am<<"   RMS_pha = "<<rms_ph<<std::endl;
			//std::cerr<<lon<<" "<<lat<<" "<<rms_ph<<"   "<<sacM.stname()<<" "<<sacS.stname()<<" "<<sacS.Dis()<<" "<<sacS.Azi()<<" rms_pha "<<sacM.fname<<"\n";
			//std::cerr<<lon<<" "<<lat<<" "<<100.*rms_am/AmpTheory<<"   "<<sacM.stname()<<" "<<sacS.stname()<<" "<<sacS.Dis()<<" "<<sacS.Azi()<<" rms_amp "<<sacM.fname<<"\n";
		}
		//std::cout<<"average misfits = "<<sqrt(amp_sum/(N-1))<<" "<<sqrt(pha_sum/(N-1))<<std::endl;
		data.Sort();

		// compute chi square
		std::vector<AziData> adVmean, adVvar;
		data.BinAverage( adVmean, adVvar, false, false );	// do not correct for 2 pi, data taken as waveform rms misfits
		for( int i=0; i<adVmean.size(); i++ ) {
			const auto& admean = adVmean[i];
			const auto& advar = adVvar[i];
			//std::cerr<<"admean = "<<admean<<"   advar = "<<advar<<"   adder = "<<(admean*admean)/advar<<std::endl;
			// (mis * mis / variance)
			chiS += adder( (admean * admean)/advar );
			N += Nadd;
		}

		// feed dataout/dataout if requested
		if( filldata ) {
			data.UpdateAziDis(minfo.lon, minfo.lat);
			dataout = std::move(data);
		}

	};

	// compute Rayleigh and Love
	chiS = 0.; N = 0;
	if( RFlag ) chiSW( _synGR, _sac3VR, dataRout );
	if( LFlag ) chiSW( _synGL, _sac3VL, dataLout );

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
			SacRec sacSZ, sacSR, sacST;
			float lon=shdM.stlo, lat=shdM.stla;
			int nptsS = ceil( (shdM.user3-minfo.t0)/shdM.delta ) + 1;
			bool synsuc = synG.ComputeSyn( sacM.stname(), lon, lat, nptsS, shdM.delta, 
													 f1,f2,f3,f4, sacSZ, sacSR, sacST, rotateSyn );
			if( ! synsuc ) {
				WarningEA::Other(FuncName, "failed to produce/invalid synthetics for station "+sacM.stname() );
				continue;
			}
			if( type == R ) {
				sacSR.clear(); sacST.clear();	// ignore unnecessory channel
			} else {
				sacSZ.clear(); sacSR.clear();
			}
			SacRec& sacS = type==R ? sacSZ : sacST;
			// check for bad sac 
			if( shdM.depmax!=shdM.depmax || sacS.shd.depmax!=sacS.shd.depmax ) continue;
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


/* -------------------- Output the data and predictions based on the input model info -------------------- */
void EQKAnalyzer::OutputFits( ModelInfo minfo ) {
	// check/correct input params
	minfo.Correct();

	// prepare/update data for the input minfo
	if( _usewaveform ) {	// fill data vectors for the waveform fitting method
		FilldataW(minfo);
	} else {	// update both path and source predictions for all SDcontainers
		PredictAll( minfo, true );
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
		sdc.BinAverage( adVmean, adVvar, isFTAN, isFTAN );	// 2pi gets corrected here
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
			const auto adstd = sqrt( adVvar[i] );
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
		FilldataW(minfo);
	} else {	// update both path and source predictions for all SDcontainers
		PredictAll( minfo, true );
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
		sdc.BinAverage( adVmean, adVvar, isFTAN, isFTAN );
		// output misfits and source preds at each bin azi (append at the end)
		AziData misL1{0.}, misL2{0.}, chiS{0.}; 
		float nbin = adVmean.size();
		for( int i=0; i<nbin; i++ ) {
			const auto& admean = adVmean[i];
			const auto& advar = adVvar[i];
			misL1 += fabs(admean); misL2 += admean * admean;
			chiS += (admean*admean) / advar;
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
	// call PredictAll for Afactors (needs a new, less redundant method here!)
	PredictAll( mi, true );
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

