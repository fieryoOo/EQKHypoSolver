#ifndef EQKANALYZER_H
#define EQKANALYZER_H

#include "DataTypes.h"
#include "SDContainer.h"
#include "FileName.h"
#include <vector>
#include <map>

class EQKAnalyzer {
public:
   /* con/destructors and operators */
   EQKAnalyzer();
   EQKAnalyzer( const std::string, bool MoveExistF = true );
   //EQKAnalyzer( const EQKAnalyzer& );
   //EQKAnalyzer( EQKAnalyzer&& );
   //EQKAnalyzer& operator= ( const EQKAnalyzer& );
   //EQKAnalyzer& operator= ( EQKAnalyzer&& );
   //~EQKAnalyzer();

   /* prepare database */
   void LoadParams( const FileName&, const bool MoveExistF = true );
   int Set( const char*, const bool MoveExistF = true );
   void CheckParams();
   void LoadData();

private:
	/* store measurements/predictions of each station with a StaData,
	 * store all StaDatas at a single period in a vector,
	 * and associate each vector with its period with a map */
	//std::map< float, std::vector<StaData> > dataR, dataL;
	std::vector<SDContainer> dataR, dataL;
};

#endif
