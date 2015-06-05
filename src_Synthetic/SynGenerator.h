#ifndef SYNGENERATOR_H
#define SYNGENERATOR_H

#include "SacRec.h"
#include "ModelInfo.h"
#include <string>

class SynGenerator {
public:
	SynGenerator( const std::string& fmodel, const ModelInfo& minfo ) 
		: fmodel(fmodel), minfo(minfo) {}
	void Run();
	bool Synthetic( const float lon, const float lat, const std::string& chname,
						 const float f1, const float f2, const float f3, const float f4, SacRec& sac );

private:
	std::string fmodel;
	ModelInfo minfo;
};


#endif
