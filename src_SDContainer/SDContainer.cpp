#include "SDContainer.h"
#include "DisAzi.h"
#include "StaList.h"
#include "VectorOperations.h"
#include <algorithm>
#include <limits>

/* IO */
void SDContainer::LoadMeasurements( const std::string& fname, const std::string fsta ) {
	// check input file
	std::ifstream fin( fname );
	if( ! fin ) 
		throw ErrorSC::BadFile(FuncName, fname);
	// clear data vector
	dataV.clear();
	// correct the input phase traveltime measurements for phase shift!
	float ph_shift = pio4_R==0 ? 0 : (-pio4_R*0.125*per);
	// read from file
	if( fsta.empty() ) {
		for(std::string line; std::getline(fin, line); ) {
			StaData sdcur(line.c_str(), ph_shift);
			if( ! sdcur.isComplete() ) continue;
			dataV.push_back( sdcur );
		}
	} else {
		StaList stalst( fsta );
		for(std::string line; std::getline(fin, line); ) {
			StaData sdcur(line.c_str(), ph_shift);
			if( ! sdcur.isComplete() ) continue;
			StaInfo data_find;
			if( stalst.SearchLoc( sdcur.lon, sdcur.lat, data_find ) )
				dataV.push_back( sdcur );
		}
	}
}

void SDContainer::LoadMaps( const std::string& fmapG, const std::string& fmapP ) {
   // read in vel maps
	mapG.Load( fmapG );
	mapP.Load( fmapP );
	if( dataV.empty() ) return;
	// get data region
	float lonmin, lonmax, latmin, latmax;
	lonmin = lonmax = dataV.at(0).lon;
	latmin = latmax = dataV.at(0).lat;
	for( const auto& sd : dataV ) {
		if( sd.lon < lonmin ) lonmin = sd.lon;
		else if( sd.lon > lonmax ) lonmax = sd.lon;
		if( sd.lat < latmin ) latmin = sd.lat;
		else if( sd.lat > latmax ) latmax = sd.lat;
	}
	std::cout<<"### Data Region: "<<lonmin<<" - "<<lonmax<<"   "<<latmin<<" - "<<latmax<<" ###\n";
	// clip maps ( for faster predictions )
	mapG.Clip( lonmin, lonmax, latmin, latmax );
	mapP.Clip( lonmin, lonmax, latmin, latmax );
}


// compute azimuth and distance for each station based on the input epicenter location
void SDContainer::UpdateAziDis( const float srclon, const float srclat ) {
	//if( lon==srclon && lat==srclat )	return;
	//lon = srclon; lat=srclat;	// do not save source location unless UpdatePathPred is called!
	for( auto& sdtmp : dataV ) {
		try {
			Path<double> pathcur(srclon, srclat, sdtmp.lon, sdtmp.lat);
			sdtmp.dis = pathcur.Dist();
			sdtmp.azi = pathcur.Azi1();
		} catch( std::exception& e ) {
			throw ErrorSC::InternalException(FuncName, e.what());
		}
		//sdtmp.valid = true;
		//if( azi<210 && azi>160 ) sdtmp.valid = false;
	}
	std::sort( dataV.begin(), dataV.end() );
}

// predict traveltimes from VelMaps and store into Gpath&Ppath
bool SDContainer::UpdatePathPred( const float srclon, const float srclat, const float srct0 ) {
	// return false if epicenter doesn't change
	if( lon==srclon && lat==srclat && t0==srct0 ) return false;
	// save new source location and update azimuth/distance
	lon = srclon; lat = srclat; t0 = srct0;
	UpdateAziDis( srclon, srclat );

	// period and wavelength
	float pero2 = per*0.5, lambda = per * Lfactor;

	if( _velG == NaN ) {
		// reset Group map and update Tpath predictions
		mapG.SetSource( srclon, srclat );
		for( auto& sd : dataV ) {
			float perc, dis = sd.dis;
			// loosen the acceptance criterion at small distances
			float lam = dis<300. ? (lambda * (2.5-dis*0.005)) : lambda;
			float minP = dis>15. ? (Min_Perc - 9./dis) : (Min_Perc - 0.6);
			//float minP = Min_Perc;
			float vel = mapG.PathAverage_Reci( Point<float>(sd.lon,sd.lat), lam, perc ).Data();
			if( perc > minP ) sd.Gpath = dis / vel + srct0;
		}
	} else {
		// update Tpaths from the fixed velocity
		for( auto& sd : dataV )	sd.Gpath = sd.dis / _velG + srct0;
	}

	if( _velP == NaN ) {
		// reset Phase map and update Tpath predictions
		mapP.SetSource( srclon, srclat );
		for( auto& sd : dataV ) {
			float perc, dis = sd.dis;
			// loosen the acceptance criterion at small distances
			float lam = dis<300. ? (lambda * (2.5-dis*0.005)) : lambda;
			float minP = dis>15. ? (Min_Perc - 9./dis) : (Min_Perc - 0.6);
			//float minP = Min_Perc;
			float vel = mapP.PathAverage_Reci( Point<float>(sd.lon,sd.lat), lam, perc ).Data();
			if( perc > minP ) sd.Ppath = dis / vel + srct0;
		}
	} else {
		for( auto& sd : dataV ) sd.Ppath = sd.dis / _velP + srct0;
	}

	return true;	// updated!
}

void SDContainer::UpdateSourcePred( const RadPattern& rad ) {
	// save new focal info
	stk = rad.stk; dip = rad.dip;
	rak = rad.rak; dep = rad.dep;
	/*
	// and update source terms
	radR.Predict( 'R', "TEST/SourceModels/245.25_41.25.R", "TEST/SourceModels/245.25_41.25.R.phv", stk, dip, rak, dep, SDContainer::perlst );
	radL.Predict( 'L', "TEST/SourceModels/245.25_41.25.L", "TEST/SourceModels/245.25_41.25.L.phv", stk, dip, rak, dep, SDContainer::perlst );
	*/

	// update source predictions
	for( auto& sd : dataV ) {
		float grt, pht, amp;
		rad.GetPred( per, sd.azi, sd.Gsource, sd.Psource, sd.Asource );
	}
}

// output a vector of AziData
//void SDContainer::ToAziVector( std::vector<AziData>& adV ) {
void SDContainer::ToMisfitV( std::vector<AziData>& adV, const float Tmin ) const {
	adV.clear();
	//const float Tmin = nwavelength * per;	// three wavelength criterion
	for( const auto& sdcur : dataV ) {
		if( sdcur.azi!=NaN ) {
			if( sdcur.Pdata < Tmin ) continue;
			AziData ad;
			if( sdcur.ToMisfit( ad ) ) adV.push_back(ad);
				//std::cerr<<"ToMisfitV:   "<<sdcur<<"   "<<ad<<std::endl;
		}
	}
}

/* compute bin average */
void SDContainer::BinAverage_ExcludeBad( std::vector<StaData>& sdVgood, bool c2pi ) {
	// dump into AziData vector
	if( c2pi ) Correct2PI();
	std::vector<AziData> adVori;
	const float Tmin = nwavelength * per;	// three wavelength criterion
	ToMisfitV( adVori, Tmin );

	// periodic extension
	std::vector<AziData> adVext;
	VO::PeriodicExtension( adVori, BINHWIDTH*2., adVext );
	adVori.clear();

	// bin average (L1 norm, [1] output with center azimuths and [2] store std-devs)
	std::vector<AziData> adVmean, adVstd;
	VO::BinAvg( adVext, adVmean, adVstd, BINSTEP, BINHWIDTH, 0., MIN_BAZI_SIZE, 1, false, false );

	// exclude empty bins and define undefined std-devs
   float finf = std::numeric_limits<float>::infinity();
	AziData ad_stdest{ BINHWIDTH/exfactor, finf, finf, finf };
	HandleBadBins(adVmean, adVstd, ad_stdest);

	// pick out good stations that falls within the range defined by adVmean and adVstd
	VO::SelectData( dataV, sdVgood, adVmean, adVstd, exfactor );
}

void SDContainer::BinAverage( std::vector<AziData>& adVmean, std::vector<AziData>& adVvar, bool c2pi, bool isFTAN ) {
	if( c2pi ) Correct2PI();
	// dump into AziData vector
	float Tmin;
	if( isFTAN ) Tmin = nwavelength * per;	// three wavelength criterion
	else Tmin = -99999.;
	std::vector<AziData> adVori;
	ToMisfitV( adVori, Tmin );
	//for( const auto& ad : adVori )	std::cerr<<ad<<std::endl;

	// periodic extension
	std::vector<AziData> adVext;
	VO::PeriodicExtension( adVori, BINHWIDTH*2., adVext );
	adVori.clear();

	// bin average (L1 norm, output with center azimuth for each bin)
	adVmean.clear(); adVvar.clear();	// note: std-devs are stored in adVvar
	VO::BinAvg( adVext, adVmean, adVvar, BINSTEP, BINHWIDTH, 0., MIN_BAZI_SIZE, 1, false, false );

	// exclude empty bins and define undefined std-devs
   float finf = std::numeric_limits<float>::infinity();
	AziData ad_stdest{ BINHWIDTH/exfactor, finf, finf, finf };
	HandleBadBins(adVmean, adVvar, ad_stdest);

	// pick out good stations that falls within the range defined by adVmean and adVvar (in which are std-devs)
	std::vector<AziData> adVsel;
	VO::SelectData( adVext, adVsel, adVmean, adVvar, exfactor );

	// bin average again (L2 norm, [1] output with center azimuths and [2] store std-dev. XX not of the means)
	VO::BinAvg( adVsel, adVmean, adVvar, BINSTEP, BINHWIDTH, 0., MIN_BAZI_SIZE, 2, true, false );

	// exclude empty bins and define undefined std-devs
	float stdest_phase = isFTAN ? stdPest : stdPHest;
	ad_stdest = AziData{ NaN, stdGest, stdest_phase, stdAest };
	// NOTE!: invalid adVvar[].Adata will be set to admean.user * ad_stdest.Adata
	HandleBadBins( adVmean, adVvar, ad_stdest );

	// compute variance by combining the data std-dev with the internally defined variance (for path predictions)
	// (old: pull up any std-devs that are smaller than defined minimum)
	AziData ad_varpath;
	if( type == R ) {
		float varmin_phase = isFTAN ? varRPmin : varRPHmin;
		ad_varpath = AziData{ NaN, varRGmin, varmin_phase, varRAmin };
	} else {
		float varmin_phase = isFTAN ? varLPmin : varLPHmin;
		ad_varpath = AziData{ NaN, varLGmin, varmin_phase, varLAmin };
	}
	//WaterLevel( adVmean, adVstd, ad_stdmin );
	ComputeVariance( adVmean, adVvar, ad_varpath );	// for adVvar: std-dev -> variance
}

// exclude empty bins and define undefined std-devs
void SDContainer::HandleBadBins( std::vector<AziData>& adVmean, 
											std::vector<AziData>& adVstd, const AziData ad_stdest ) const {
   std::vector<AziData> adVmean2, adVstd2;
   for(int i=0; i<adVmean.size(); i++) {
      auto& admean = adVmean[i];
      auto& adstd = adVstd[i];
      // continue if the current bin was empty
      if( admean.azi == NaN ) continue;
      // re-assign azi range and data std-devs
      if( ad_stdest.azi != NaN ) adstd.azi = ad_stdest.azi;
      if( adstd.Gdata == NaN ) adstd.Gdata = ad_stdest.Gdata;
      if( adstd.Pdata == NaN ) adstd.Pdata = ad_stdest.Pdata;
      if( adstd.Adata == NaN ) adstd.Adata = admean.user * ad_stdest.Adata;
      // store
      adVmean2.push_back( std::move(admean) );
      adVstd2.push_back( std::move(adstd) );
      //std::cerr<<"mean = "<<admean<<"   std = "<<adstd<<std::endl;
   }
	adVmean = std::move( adVmean2 );
	adVstd  = std::move( adVstd2  );
}

// compute variance by propagating the given variance into the data
// ( old: pull up any std-devs that are smaller than defined minimum )
void SDContainer::ComputeVariance( std::vector<AziData>& adVmean, 
											  std::vector<AziData>& adVstd, const AziData ad_varin ) const {
	float Gin = ad_varin.Gdata, Pin = ad_varin.Pdata, Asca = ad_varin.Adata;
   for(int i=0; i<adVmean.size(); i++) {
      auto& admean = adVmean[i];
      auto& adstd = adVstd[i];
//std::cerr<<"before var correction: admean = "<<admean<<"   adstd = "<<adstd<<"   admean.user = "<<admean.user<<" Aperc = "<<Asca<<std::endl;
		//if( adstd.Gdata < Gmin ) adstd.Gdata = Gmin;
		//if( adstd.Pdata < Pmin ) adstd.Pdata = Pmin;
		//float Amin = admean.user * Asca;
		//if( adstd.Adata < Amin ) adstd.Adata = Amin;
		float Ain = admean.user*admean.user * Asca;
		adstd.Gdata = adstd.Gdata*adstd.Gdata + Gin;
		adstd.Pdata = adstd.Pdata*adstd.Pdata + Pin;
		adstd.Adata = adstd.Adata*adstd.Adata + Ain;
//std::cerr<<"after var correction: admean = "<<admean<<"   advar = "<<adstd<<std::endl;
	}
}

