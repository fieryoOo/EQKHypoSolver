#include "DataTypes.h"
#include "EQKAnalyzer.h"
#include "VectorOperations.h"
#include "SDContainer.h"
#include <iostream>
#include <fstream>

int main( int argc, char* argv[] ) {
	if( argc != 4 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [infile] [mapG] [mapP]"<<std::endl;
		return -1;
	}

	// load data
	float per = 22.;
	SDContainer sdc( per, argv[1], argv[2], argv[3] );

	const float srclon = -114.8, srclat = 41.26;

	// predict travel time along GCPs
	sdc.UpdatePathPred( srclon, srclat );

	// update source term for each station
	const float stk = 210, dip = 40, rak = -90., dep = 8.;
	RadPattern radR, radL;
	std::vector<float> perlst{ per };
	radR.Predict( 'R', "TEST/SourceModels/245.25_41.25.R", "TEST/SourceModels/245.25_41.25.R.phv", stk, dip, rak, dep, perlst );
	radL.Predict( 'L', "TEST/SourceModels/245.25_41.25.L", "TEST/SourceModels/245.25_41.25.L.phv", stk, dip, rak, dep, perlst );

	sdc.UpdateSourcePred( radR );
	std::ofstream fout("dat.sta");
	for( const auto& sd : sdc.dataV )
		fout<<sd<<"\n";
	fout.close(); fout.clear();

	// scale source amplitudes to match the observations
	float Afactor = 50000.;
	sdc.AmplifySource( Afactor );

	// compute misfits on all stations,
	// average in each (20 degree) bin,
	// and store the results into AziData vectors
	std::vector<AziData> adVmean, adVstd;
	sdc.BinAverage( adVmean, adVstd );
	fout.open("res.azi");
	for(int i=0; i<adVmean.size(); i++) {
		const auto& admean = adVmean[i];
		const auto& adstd = adVstd[i];
		fout<<admean<<" "<<adstd<<std::endl;
	}
	fout.close(); fout.clear();

	std::vector<StaData> sdVgood;
	sdc.BinAverage_ExcludeBad( sdVgood );
	fout.open("res.sta");
	for( const auto& sd : sdVgood )
		fout<<sd<<"\n";

	return 0;
}
