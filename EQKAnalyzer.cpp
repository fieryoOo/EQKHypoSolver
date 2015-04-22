#include "EQKAnalyzer.h"
#include "DisAzi.h"
#include <sstream>
#include <sys/stat.h>

int WarningEA::Base::nWarnOther = 0;

EQKAnalyzer::EQKAnalyzer() {}

EQKAnalyzer::EQKAnalyzer( const std::string fparam, bool MoveF ) {
	LoadParams( fparam, MoveF );
	CheckParams();
}

//EQKAnalyzer::~EQKAnalyzer() {}


/* -------------------- param/data preparations (loading) -------------------- */
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
		_model.lon = clon;
		_model.lat = clat;
		_model.t0 = ct0;
		return true;
	}
	else {
		Rlon = Rlat = NaN;
		return false;
	}
}

bool EQKAnalyzer::MKDir(const char *dirname) {
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

void EQKAnalyzer::MKDirFor( const std::string& path ) {
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

   std::cout<<"### "<<nparam<<" succed loads from param file "<<fname<<". ###"<<std::endl;

   // check clon/clat validation and initialize einfo
   if( ! InitEpic() ) throw ErrorEA::BadParam(FuncName, "invalid Rs | clon | clat");

   // synchronize output vectors with input-measurement vectors
   //pimplES->ArrangeOutlist( pimplES->outlist_RF, pimplES->perRtmp_fit, pimplES->perRlist );
   //pimplES->ArrangeOutlist( pimplES->outlist_LF, pimplES->perLtmp_fit, pimplES->perLlist );
   //pimplES->ArrangeOutlist( pimplES->outlist_RP, pimplES->perRtmp_pred, pimplES->perRlist );
   //pimplES->ArrangeOutlist( pimplES->outlist_LP, pimplES->perLtmp_pred, pimplES->perLlist );

   // check output paths and make directorie(s) if necessary
   std::vector<FileName> Voutf { outname_misF, outname_misL, outname_misAll, outname_pos };
	for( auto fname : outlist_RF ) Voutf.push_back( fname.second );
	for( auto fname : outlist_LF ) Voutf.push_back( fname.second );
	for( auto fname : outlist_RP ) Voutf.push_back( fname.second );
	for( auto fname : outlist_LP ) Voutf.push_back( fname.second );
   for( const auto& outfname : Voutf ) MKDirFor( outfname );
}

int EQKAnalyzer::Set( const char *input, const bool MoveF ) {
	std::istringstream buff(input);
	std::string stmp;
	if( ! (buff>>stmp) ) return -2;
	bool succeed;
	if( stmp == "clon" ) { 
		succeed = buff >> clon; 
		if( succeed && clon<0.) clon += 360.; 
	}
	else if( stmp == "clat" ) succeed = buff >> clat;
	else if( stmp == "ct0" ) succeed = buff >> ct0;
	else if( stmp == "Rs") succeed = buff >> Rs;
	else if( stmp == "strike" ) succeed = buff >> _model.strike;
	else if( stmp == "dip" ) succeed = buff >> _model.dip;
	else if( stmp == "rake" ) succeed = buff >> _model.rake;
	else if( stmp == "depth" ) succeed = buff >> _model.depth;
	else if( stmp == "fRse" ) succeed = buff >> fReigname;
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
		if( succeed && MoveF && access(outname.c_str(), F_OK) == 0 ) {
			FileName oldname = outname + "_old";
			WarningEA::MoveExistFile(FuncName, outname+" -> "+oldname);
			rename(outname.c_str(), oldname.c_str());
		}
	}
	else if( stmp == "fmisF" ) {
		FileName& outname = outname_misF;
		succeed = buff >> outname;
		if( succeed && MoveF && access(outname.c_str(), F_OK) == 0 ) {
			FileName oldname = outname + "_old";
			WarningEA::MoveExistFile(FuncName, outname+" -> "+oldname);
			rename(outname.c_str(), oldname.c_str());
		}
	}
	else if( stmp == "fmisAll" ) {
		FileName& outname = outname_misAll;
		succeed = buff >> outname;
		if( succeed && MoveF && access(outname.c_str(), F_OK) == 0 ) {
			FileName oldname = outname + "_old";
			WarningEA::MoveExistFile(FuncName, outname+" -> "+oldname);
			rename(outname.c_str(), oldname.c_str());
		}
	}
	else if( stmp == "fpos" ) {
		FileName& outname = outname_pos;
		succeed = buff >> outname;
		if( succeed && MoveF && access(outname.c_str(), F_OK) == 0 ) {
			FileName oldname = outname + "_old";
			WarningEA::MoveExistFile(FuncName, outname+" -> "+oldname);
			rename(outname.c_str(), oldname.c_str());
		}
	}
	else if( stmp == "ffitR" ) {
		FileName outname;
		float per;
		succeed = buff >> outname >> per;
		if( succeed ) {
			if( MoveF && access(outname.c_str(), F_OK) == 0 ) {
				FileName oldname = outname + "_old";
				WarningEA::MoveExistFile(FuncName, outname+" -> "+oldname);
				rename(outname.c_str(), oldname.c_str());
			}
			outlist_RF[per] = outname;
		}
	}
	else if( stmp == "ffitL" ) {
		FileName outname;
		float per;
		succeed = buff >> outname >> per;
		if( succeed ) {
			if( MoveF && access(outname.c_str(), F_OK) == 0 ) {
				FileName oldname = outname + "_old";
				WarningEA::MoveExistFile(FuncName, outname+" -> "+oldname);
				rename(outname.c_str(), oldname.c_str());
			}
			outlist_LF[per] = outname;
		}
	}
	else if( stmp == "fpredR" ) {
		FileName outname;
		float per;
		succeed = buff >> outname >> per;
		if( succeed ) {
			/*
			if( MoveF && access(outname.c_str(), F_OK) == 0 ) {
				FileName oldname = outname + "_old";
				WarningEA::MoveExistFile(FuncName, outname+" -> "+oldname);
				rename(outname.c_str(), oldname.c_str());
			} */
			outlist_RP[per] = outname;
		}
	}
	else if( stmp == "fpredL" ) {
		FileName outname;
		float per;
		succeed = buff >> outname >> per;
		if( succeed ) {
			/*
			if( MoveF && access(outname.c_str(), F_OK) == 0 ) {
				FileName oldname = outname + "_old";
				WarningEA::MoveExistFile(FuncName, outname+" -> "+oldname);
				rename(outname.c_str(), oldname.c_str());
			} */
			outlist_LP[per] = outname;
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
	const float NaN = ModelInfo::NaN;
	if( _model.strike==NaN|| _model.dip==NaN||
		 _model.rake==NaN|| _model.depth==NaN) {
      throw ErrorEA::BadParam(FuncName, "invalid/incomplete initial focal info");
   }

   // check model inputs
   if( access(fReigname.c_str(), F_OK) == -1 ) throw ErrorEA::BadFile(FuncName, fReigname);
   if( access(fRphvname.c_str(), F_OK) == -1 ) throw ErrorEA::BadFile(FuncName, fRphvname);
   if( access(fLeigname.c_str(), F_OK) == -1 ) throw ErrorEA::BadFile(FuncName, fLeigname);
   if( access(fLphvname.c_str(), F_OK) == -1 ) throw ErrorEA::BadFile(FuncName, fLphvname);

	// normalize data weightings
	NormalizeWeights( datatype, weightR_Loc, weightL_Loc );
	NormalizeWeights( datatype, weightR_Foc, weightL_Foc );

	// check excluded data
   if( !_useG && !_useP )
      throw ErrorEA::BadParam(FuncName, "All data excluded! (noG and noP)");
   if( ! _useG )
      WarningEA::Other(FuncName, "GroupT data excluded!" );
   if( ! _useP )
      WarningEA::Other(FuncName, "PhaseT data excluded!" );
   if( ! _useA )
      WarningEA::Other(FuncName, "Amplitude data excluded!" );

}


/* -------------------- Check/Load all data file requested by fparam -------------------- */
void EQKAnalyzer::LoadData() {
	// Rayleigh
	_dataR.clear();
	for( const auto& fR : fRlist ) {
		const float per = fR.first;
		const auto& farray = fR.second;
		_dataR.push_back( SDContainer(per, farray[0], farray[1], farray[2]) );
	}

	// Love
   _dataL.clear(); 
	for( const auto& fL : fLlist ) {
		const float per = fL.first;
		const auto& farray = fL.second;
		_dataL.push_back( SDContainer(per, farray[0], farray[1], farray[2]) );
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

void EQKAnalyzer::InitRadPattern() {
	float stk = _model.strike, dip = _model.dip, rak = _model.rake, dep = _model.depth;
	_rpR.Predict( 'R', fReigname, fRphvname, stk, dip, rak, dep, perRlst() );
	_rpL.Predict( 'L', fLeigname, fLphvname, stk, dip, rak, dep, perLlst() );
}

/* -------------------- compute the total chi-square misfit based on the current data state and the input model info -------------------- */
void EQKAnalyzer::chiSquare( const ModelInfo& minfo, float& chiS, float& wSum, int& N, bool isInit ) const {
	// check input params
	if( ! minfo.isValid() )
		throw ErrorEA::BadParam(FuncName, "invalid model parameter(s)");
	int Rsize = _dataR.size(), Lsize = _dataL.size();
	if( Rsize==0 && Lsize==0 )
		throw ErrorEA::EmptyData(FuncName, "Rsize && Lsize");
	if( datatype == Undefined )
		throw ErrorEA::BadParam(FuncName, "Undefined datatype");

	// data flags
	bool useG = _useG, useP = _useP, useA = _useA;
	bool RFlag = (datatype==B || datatype==R);
	bool LFlag = (datatype==B || datatype==L);
	if( isInit ) {
		useG = LFlag = true;
		useP = RFlag = false;
	}

	// model info
	float stk = minfo.strike, dip = minfo.dip, rak = minfo.rake, dep = minfo.depth;

	// main loop for Rayleigh
	chiS = 0.; wSum = 0.; N = 0;
	if( Rsize>0 && RFlag ) {
		// update source term predictions
		auto rpR = _rpR;
		rpR.Predict( 'R', fReigname, fRphvname, stk, dip, rak, dep, perRlst() );
		// make a copy at each period (for multi-threading)
		for( auto sdc : _dataR ) {
			sdc.UpdatePathPred( minfo.lon, minfo.lat );
			sdc.UpdateSourcePred( rpR );
			float Afactor = 50000.;
			// scale source amplitudes to match the observations
			sdc.AmplifySource( Afactor );
			// compute misfits on all stations,
			// average in each (20 degree) bin,
			// and store the results into AziData vectors
			std::vector<AziData> adVmean, adVstd;
			sdc.BinAverage( adVmean, adVstd );
for( const auto& ad : adVmean )
	std::cerr<<ad<<"\n";
for( const auto& ad : adVstd )
	std::cerr<<ad<<"\n";
std::cerr<<"chiSquare: per = "<<sdc.per<<" done!\n";
			//std::vector<StaData> sdVgood;
			//sdc.BinAverage_ExcludeBad( sdVgood );
			// compute misfit
		}
	}

}
