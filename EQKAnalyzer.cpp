#include "EQKAnalyzer.h"
#include "DisAzi.h"
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
/*
bool EQKAnalyzer::InitEpic() {
	const float NaN = ModelInfo::NaN;
	if( Rs>0 && clon!=NaN && clat!=NaN) {
		double dist;
		try {
			//calc_dist(clat, clon, clat, clon+1, &dist);
			dist = Path<double>(clon, clat, clon+1., clat).Dist();
			Rlon = Rs / dist;
			//calc_dist(clat, clon, clat+1, clon, &dist);
			dist = Path<double>(clon, clat, clon, clat+1.).Dist();
			Rlat = Rs / dist;
		} catch ( std::exception& e ) {
			std::cerr<<e.what()<<std::endl;
			return false;
		}
		return true;
	}
	else {
		Rlon = Rlat = NaN;
		return false;
	}
}
*/

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

void EQKAnalyzer::MKDirFor( const std::string& path ) const {
	if( path.empty() ) return;
	std::stringstream spath(path);
	std::vector<std::string> dirs;
	for(std::string dir; std::getline(spath, dir, '/'); ) {
		dirs.push_back(dir);
	}
	if( dirs.size() < 2 ) return;
	for( int idir=0; idir<dirs.size()-1; idir++ ) {
		if( MKDir( dirs[idir].c_str() ) )
			WarningEA::Other(FuncName, "Making directory " + dirs[idir] );
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

   std::cout<<"### EQKAnalyzer::LoadParams: "<<nparam<<" succed loads from param file "<<fname<<". ###"<<std::endl;

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
	for( auto fname : outlist_RP ) Voutf.push_back( fname.second );
	for( auto fname : outlist_LP ) Voutf.push_back( fname.second );
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
		FileName stmp1, stmp2, stmp3;
		float per;
		succeed = buff >> stmp1 >> stmp2 >> stmp3 >> per;
		if( succeed ) {
			if( fRlist.find(per) != fRlist.end() ) return -3;	// period already exists
			fRlist[per] = std::array<FileName, 3>{stmp1, stmp2, stmp3};
		}
	}
	else if( stmp == "fLm" ) {
		FileName stmp1, stmp2, stmp3;
		float per;
		succeed = buff >> stmp1 >> stmp2 >> stmp3 >> per;
		if( succeed ) {
			if( fLlist.find(per) != fLlist.end() ) return -3;	// period already exists
			fLlist[per] = std::array<FileName, 3>{stmp1, stmp2, stmp3};
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
	}
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
	// Rayleigh
	_dataR.clear();
	for( const auto& fR : fRlist ) {
		const float per = fR.first;
		const auto& farray = fR.second;
		// farray[0]: measurement file
		// farray[1] & [2]: either G&P vel_map files or G&P velocities (float for 1D model)
		float velG, velP;
		if( FilenameToVel(farray[1], velG) && FilenameToVel(farray[2], velP) ) {
			_dataR.push_back( SDContainer(per, 'R', farray[0], velG, velP) );
		} else {
			_dataR.push_back( SDContainer(per, 'R', farray[0], farray[1], farray[2]) );
		}
	}

	// Love
   _dataL.clear(); 
	for( const auto& fL : fLlist ) {
		const float per = fL.first;
		const auto& farray = fL.second;
		float velG, velP;
		if( FilenameToVel(farray[1], velG) && FilenameToVel(farray[2], velP) ) {
			_dataL.push_back( SDContainer(per, 'L', farray[0], velG, velP) );
		} else {
			_dataL.push_back( SDContainer(per, 'L', farray[0], farray[1], farray[2]) );
		}
	}

   std::cout<<"### "<<_dataR.size() + _dataL.size()<<" data file(s) loaded. ###"<<std::endl;
}


/* -------------------- fill the rpR and rpL objects with the current model state -------------------- */
std::vector<float> EQKAnalyzer::perRlst() const {
	std::vector<float> perlst;
	for( const auto& data : _dataR ) perlst.push_back( data.per );
	return perlst;
}
std::vector<float> EQKAnalyzer::perLlst() const {
	std::vector<float> perlst;
	for( const auto& data : _dataL ) perlst.push_back( data.per );
	return perlst;
}

// initialize the Analyzer by pre- predicting radpatterns and updating pathpred for all SDContainer based on the current model info
void EQKAnalyzer::PredictAll( const ModelInfo& mi, bool updateSource ) {
	PredictAll( mi, _rpR, _rpL, _dataR, _dataL, _AfactorR, _AfactorL, _source_updated, updateSource );
}
void EQKAnalyzer::PredictAll( const ModelInfo& mi,	std::vector<SDContainer>& dataR, 
										std::vector<SDContainer>& dataL, bool updateSource ) const {
	auto rpR = _rpR, rpL = _rpL;
	float AfactorR = _AfactorR, AfactorL = _AfactorL;
	bool source_updated = _source_updated;
	PredictAll( mi, rpR, rpL, dataR, dataL, AfactorR, AfactorL, source_updated, updateSource );
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
void EQKAnalyzer::PredictAll( const ModelInfo& mi, RadPattern& rpR, RadPattern& rpL,
										std::vector<SDContainer>& dataR, std::vector<SDContainer>& dataL,
										float& AfactorR, float& AfactorL, bool& source_updated, bool updateSource ) const {
	// radpattern
	bool model_updated = false;
	float stk = mi.stk, dip = mi.dip, rak = mi.rak, dep = mi.dep;
	stk = ShiftInto( stk, 0., 360., 360. );	dip = BoundInto( dip, 0., 90. );
	rak = ShiftInto( rak, -180., 180., 360. ); dep = BoundInto( dep, 0., 60. );
	model_updated |= rpR.Predict( 'R', fReigname, fRphvname, stk, dip, rak, dep, perRlst() );
	model_updated |= rpL.Predict( 'L', fLeigname, fLphvname, stk, dip, rak, dep, perLlst() );

	// SDContainer Tpaths
	for( auto& sdc : dataR )
		model_updated |= sdc.UpdatePathPred( mi.lon, mi.lat, mi.t0 );
	for( auto& sdc : dataL )
		model_updated |= sdc.UpdatePathPred( mi.lon, mi.lat, mi.t0 );

	if( model_updated ) source_updated = false;	// source terms need to be updated
	else if( source_updated ) return;				//	model state didn't change since the last source term update -- return
	if( ! updateSource ) return;						// do not update source for now
	source_updated = true;								// do update source (now)

	// SDContainer R: source terms
	if( dataR.size() > 0 ) {
		std::vector<float> ampratioV;	
		ampratioV.reserve( dataR.size() * dataR.at(0).size() );
		for( auto& sdc : dataR ) {
			sdc.UpdateSourcePred( rpR );
			sdc.ComputeAmpRatios( ampratioV );
		}
		AfactorR = std::accumulate(ampratioV.begin(), ampratioV.end(), 0.) / ampratioV.size();
		//std::cerr<<"   Afactor for Rayleigh: "<<Afactor<<" (was 50000.)\n";
		for( auto& sdc : dataR )	// scale source amplitudes to match the observations
			sdc.AmplifySource( AfactorR );
	}

	// SDContainers L: source terms
	if( dataL.size() > 0 ) {
		std::vector<float> ampratioV;	
		ampratioV.reserve( dataL.size() * dataL.at(0).size() );
		for( auto& sdc : dataL ) {
			sdc.UpdateSourcePred( rpL );
			sdc.ComputeAmpRatios( ampratioV );
		}
		AfactorL = std::accumulate(ampratioV.begin(), ampratioV.end(), 0.) / ampratioV.size();
		//std::cerr<<"   Afactor for Love: "<<Afactor<<" (was 55000.)\n";
		for( auto& sdc : dataL )	// scale source amplitudes to match the observations
			sdc.AmplifySource( AfactorL );
	}
}

/* -------------------- compute the total chi-square misfit based on the current data state and the input model info -------------------- */
void EQKAnalyzer::chiSquare( const ModelInfo& minfo, float& chiS, int& N ) const {
	// check input params
	if( ! minfo.isValid() ) {
		std::stringstream ss; ss<<minfo;
		throw ErrorEA::BadParam( FuncName, "invalid model parameter(s): " + ss.str() );
	}
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
		useG = LFlag = true;
		useA = useP = RFlag = false;
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
			sdc.BinAverage( adVmean, adVvar );
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
			sdc.BinAverage( adVmean, adVvar );
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


/* -------------------- Output the data and predictions based on the input model info -------------------- */
void EQKAnalyzer::Output( const ModelInfo& minfo ) {
	// check input params
	if( ! minfo.isValid() ) {
		std::stringstream ss; ss<<minfo;
		throw ErrorEA::BadParam( FuncName, "invalid model parameter(s): " + ss.str() );
	}
	int Rsize = _dataR.size(), Lsize = _dataL.size();
	if( Rsize==0 && Lsize==0 )
		throw ErrorEA::EmptyData(FuncName, "Rsize && Lsize");

	// make a copy of the data and RadPattern objs (for multi-threading),
	// update both path and source predictions for all SDcontainers
	//auto rpR(_rpR), rpL(_rpL);
	//auto dataR(_dataR), dataL(_dataL);
	//bool source_updated = _source_updated;
	//PredictAll( minfo, rpR, rpL, dataR, dataL, source_updated, true );
	PredictAll( minfo, true );

	// main loop for Rayleigh
	if( Rsize>0 ) {
		for( auto& sdc : _dataR ) {
			const float per = sdc.per;
			// outname for the current period
			if( outlist_RF.find(per) == outlist_RF.end() ) continue;
			const std::string& outname = outlist_RF.at(per);
			// average in each (20 degree) bin,
			std::vector<AziData> adVmean, adVvar;
			sdc.BinAverage( adVmean, adVvar );
			// output data and predictions at each station (append at the end)
			std::ofstream foutsta( outname + "_sta", std::ofstream::app );
			foutsta<<"# [ minfo = "<<minfo<<" ]\n";
			sdc.PrintAll( foutsta );
			foutsta << "\n\n";
			// output misfits and source preds at each bin azi (append at the end)
			std::ofstream foutbin( outname + "_bin", std::ofstream::app );
			foutbin<<"# [ minfo = "<<minfo<<" ]\n";
			for( int i=0; i<adVmean.size(); i++ ) {
				const auto& admean = adVmean[i];
				const auto adstd = sqrt( adVvar[i] );
				float Gsource, Psource, Asource;
				_rpR.GetPred( per, admean.azi,	Gsource, Psource, Asource );
				Asource *= _AfactorR;
				foutbin<<admean.azi<<"  "<<Gsource<<" "<<admean.Gdata<<" "<<adstd.Gdata
										 <<"  "<<Psource<<" "<<admean.Pdata<<" "<<adstd.Pdata
										 <<"  "<<Asource<<" "<<admean.Adata<<" "<<adstd.Adata<<"\n";
			}
			foutbin << "\n\n";
		}
	}

	// main loop for Love
	if( Lsize>0 ) {
		for( auto& sdc : _dataL ) {
			const float per = sdc.per;
			// outname for the current period
			if( outlist_LF.find(per) == outlist_LF.end() ) continue;
			const std::string& outname = outlist_LF.at(per);
			// average in each (20 degree) bin,
			std::vector<AziData> adVmean, adVvar;
			sdc.BinAverage( adVmean, adVvar );
			// output data and predictions at each station (append at the end)
			std::ofstream foutsta( outname + "_sta", std::ofstream::app );
			foutsta<<"# [ minfo = "<<minfo<<" ]\n";
			sdc.PrintAll( foutsta );
			foutsta << "\n\n";
			// output misfits and source preds at each bin azi (append at the end)
			std::ofstream foutbin( outname + "_bin", std::ofstream::app );
			foutbin<<"# [ minfo = "<<minfo<<" ]\n";
			for( int i=0; i<adVmean.size(); i++ ) {
				const auto& admean = adVmean[i];
				const auto adstd = sqrt( adVvar[i] );
				float Gsource, Psource, Asource;
				_rpL.GetPred( per, admean.azi,	Gsource, Psource, Asource );
				Asource *= _AfactorL;
				foutbin<<admean.azi<<"  "<<Gsource<<" "<<admean.Gdata<<" "<<adstd.Gdata
										 <<"  "<<Psource<<" "<<admean.Pdata<<" "<<adstd.Pdata
										 <<"  "<<Asource<<" "<<admean.Adata<<" "<<adstd.Adata<<"\n";
			}
			foutbin << "\n\n";
		}
	}

}


/* -------------------- compute and output misfits, separately, for group, phase, and amplitudes -------------------- */
void EQKAnalyzer::OutputMisfits( const ModelInfo& minfo ) {
	// check input params
	if( ! minfo.isValid() ) {
		std::stringstream ss; ss<<minfo;
		throw ErrorEA::BadParam( FuncName, "invalid model parameter(s): " + ss.str() );
	}
	int Rsize = _dataR.size(), Lsize = _dataL.size();
	if( Rsize==0 && Lsize==0 )
		throw ErrorEA::EmptyData(FuncName, "Rsize && Lsize");

	PredictAll( minfo, true );

	// prepare output file
	std::ofstream fout( outname_misAll, std::ofstream::app );
	fout<<"# [ minfo = "<<minfo<<" ]\n";
	fout<<"wtype per  rmsG L1G rchisG  rmsP L1P rchisP  rmsA L1A rchisA\n";

	// main loop for Rayleigh
	if( Rsize>0 ) {
		for( auto& sdc : _dataR ) {
			const float per = sdc.per;
			// average in each (20 degree) bin,
			std::vector<AziData> adVmean, adVvar;
			sdc.BinAverage( adVmean, adVvar );
			// output misfits and source preds at each bin azi (append at the end)
			AziData misL1{0.}, misL2{0.}, chiS{0.}; 
			float nbin = adVmean.size();
			for( int i=0; i<nbin; i++ ) {
				const auto& admean = adVmean[i];
				const auto& advar = adVvar[i];
				misL1 += fabs(admean); misL2 += admean * admean;
				chiS += (admean*admean) / advar;
/*
if( per == 22. ) {
	std::cerr<<" *** debug *** for R at per="<<per<<" azi="<<admean.azi<<": misA="<<admean.Adata<<" stdA="<<sqrt(advar.Adata)<<" chiSAsum="<<chiS.Adata<<"\n";
}
*/
			}
			misL1 /= nbin; misL2 = sqrt( misL2/(nbin-1.) );
			chiS /= nbin;
			fout << "R " << per << "  "
				  << misL2.Gdata << " " << misL1.Gdata << " " <<  chiS.Gdata << "  "
				  << misL2.Pdata << " " << misL1.Pdata << " " <<  chiS.Pdata << "  "
				  << misL2.Adata << " " << misL1.Adata << " " <<  chiS.Adata << "\n";
		}
	}

	// main loop for Love
	if( Lsize>0 ) {
		for( auto& sdc : _dataL ) {
			const float per = sdc.per;
			// average in each (20 degree) bin,
			std::vector<AziData> adVmean, adVvar;
			sdc.BinAverage( adVmean, adVvar );
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
			fout << "L " << per << "  "
				  << misL2.Gdata << " " << misL1.Gdata << " " <<  chiS.Gdata << "  "
				  << misL2.Pdata << " " << misL1.Pdata << " " <<  chiS.Pdata << "  "
				  << misL2.Adata << " " << misL1.Adata << " " <<  chiS.Adata << "\n";
		}
	}

}
