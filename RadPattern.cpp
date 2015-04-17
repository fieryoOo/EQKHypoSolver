#include "RadPattern.h"
#include <fstream>

/* FORTRAN entrance */
const int nazi = RadPattern::nazi;
extern"C" {
   void rad_pattern_r_(char *feig_buff, int *eig_len, int *phvnper, int *phvdper,
                       const float *strike, const float *dip, const float *rake, const float *depth,
		       const float *per, int *nper,
                       float *azi, float grT[][nazi], float phT[][nazi], float amp[][nazi]);
   void rad_pattern_l_(char *feig_buff, int *eig_len, int *phvnper, int *phvdper,
                       const float *strike, const float *dip, const float *rake, const float *depth,
		       const float *per, int *nper,
                       float *azi, float grT[][nazi], float phT[][nazi], float amp[][nazi]);
}


/* ---------- implimentation ---------- */
struct RadPattern::Rimpl {
   //FocalInfo<ftype> finfo, finfoold;
   //std::string outname_mis;
   std::string feignmem, fphvnmem; // name of files of the currently-in-the-memory contents

   int phvnper, phvdper;
   int feig_len;
   char *feig_buff;

   /* ---------- con/destructors ---------- */
   Rimpl() {
      phvnper = -12345; phvdper = -12345;
      feig_buff = nullptr; feig_len = -12345;
   }
   ~Rimpl() { if(feig_buff) delete [] feig_buff; }

};


/* con/destructors and operators */
RadPattern::RadPattern()
   : pimplR( new Rimpl ) {}

RadPattern::RadPattern( const RadPattern& rp2 )
   : pimplR( new Rimpl(*(rp2.pimplR)) ) {}

RadPattern::RadPattern( RadPattern&& rp2 )
   : pimplR( std::move(rp2.pimplR) ) {}

RadPattern& RadPattern::operator= ( const RadPattern& rp2 ) {
   pimplR.reset( new Rimpl(*(rp2.pimplR)) );
}

RadPattern& RadPattern::operator= ( RadPattern&& rp2 ){
   pimplR = std::move(rp2.pimplR);
}

RadPattern::~RadPattern() {}

/* predict radpattern for rayleigh and love waves */
void RadPattern::Predict( char typein, const std::string& feigname, const std::string& fphvname,
			  const ftype stkin, const ftype dipin, const ftype rakin, const ftype depin,
			  const std::vector<float>& perlst ) {
			  //std::vector< std::vector<AziData> >& per_azi_pred ) {
   if( typein!='R' && typein!='L' )
      throw ErrorRP::BadParam(FuncName, "unknown type = "+type);

	// return if the requested new state is exactly the same as the one stored
	if( type==typein && stk==stkin && dip==dipin && 
		 rak==rakin && dep==depin && perlst.size()<=grtM.size() ) {
		for( const auto per : perlst )
			if( grtM.find(per) == grtM.end() ) return;
	}

	// store current state;
	type = typein; 
	stk = stkin; dip = dipin;
	rak = rakin; dep = depin;

  #pragma omp critical
  { // omp critical begins
   // read feig into memory
   if( pimplR->feig_buff == nullptr || feigname != pimplR->feignmem ) {
      if( feigname.empty() ) throw ErrorRP::BadParam(FuncName, "empty feigname");
      std::ifstream fin( feigname.c_str() );
      if( ! fin ) throw ErrorRP::BadFile(FuncName, feigname);
      fin.seekg(0, std::ios::end);
      pimplR->feig_len = fin.tellg();
      if( pimplR->feig_buff ) {
	 delete [] pimplR->feig_buff;
	 pimplR->feig_buff = nullptr;
      }
      pimplR->feig_buff = new char[pimplR->feig_len];
      fin.seekg(0,std::ios::beg);
      fin.read(pimplR->feig_buff, pimplR->feig_len);
      fin.close();
      pimplR->feignmem = feigname;
   }

   // read in nper and dper from fphv
   if( pimplR->phvnper < 0 || fphvname != pimplR->fphvnmem )  {
      if( fphvname.empty() ) throw ErrorRP::BadParam(FuncName, "empty fphvname");
      std::ifstream fin( fphvname.c_str() );
      if( ! fin ) throw ErrorRP::BadFile(FuncName, fphvname);
      pimplR->phvnper = 0;
      for( std::string line; std::getline(fin, line); ) {
         float ftmp1, ftmp2, ftmp3;
         if( sscanf(line.c_str(), "%f %f %f", &ftmp1, &ftmp2, &ftmp3) != 3 ) continue;
         if( pimplR->phvnper == 0 ) pimplR->phvdper = ftmp1;
         else if( pimplR->phvnper == 1 ) pimplR->phvdper = ftmp1 - pimplR->phvdper;
         (pimplR->phvnper)++;
      }
      fin.close();
      pimplR->fphvnmem = fphvname;
   }
  } // omp critical ends

   // check if feig and fphv contents are modified
   if( pimplR->feig_buff == nullptr || feigname != pimplR->feignmem ) 
      throw ErrorRP::BadBuff(FuncName, feigname + " != " + pimplR->feignmem);
   if( pimplR->phvnper < 0 || fphvname != pimplR->fphvnmem ) 
      throw ErrorRP::BadBuff(FuncName, fphvname + " != " + pimplR->fphvnmem);

   // run rad_pattern
   int nper = perlst.size();
	float azi[nazi], grT[nper][nazi], phT[nper][nazi], amp[nper][nazi];

   if( type == 'R' ) {
      rad_pattern_r_( pimplR->feig_buff, &(pimplR->feig_len), &(pimplR->phvnper), &(pimplR->phvdper),
                      &(stk), &(dip), &(rak), &(dep), &(perlst.at(0)), &nper, azi, grT, phT, amp );
   } else if( type == 'L' ) {
      rad_pattern_l_( pimplR->feig_buff, &(pimplR->feig_len), &(pimplR->phvnper), &(pimplR->phvdper),
                      &(stk), &(dip), &(rak), &(dep), &(perlst.at(0)), &nper, azi, grT, phT, amp );
   }

   // check if feig and fphv contents are modified
   if( pimplR->feig_buff == nullptr || feigname != pimplR->feignmem ) 
      throw ErrorRP::BadBuff(FuncName, feigname + " != " + pimplR->feignmem);
   if( pimplR->phvnper < 0 || fphvname != pimplR->fphvnmem ) 
      throw ErrorRP::BadBuff(FuncName, fphvname + " != " + pimplR->fphvnmem);

   // copy results into prediction vectors
/*
   per_azi_pred.clear(); per_azi_pred.resize(nper);
   for(int i=0; i<nper; i++) {
      per_azi_pred.at(i).resize(nazi);
      for(int iazi=0; iazi<nazi; iazi++) {
         if( azi[iazi] != iazi*2 )
	    throw ErrorRP::BadAzi(FuncName, std::to_string(azi[iazi])+" != 2*"+std::to_string(iazi));
         per_azi_pred[i][iazi] = AziData( azi[iazi], grT[i][iazi], phT[i][iazi], amp[i][iazi] );
      }
   }
*/
   //float strike = finfo.strike, dip = finfo.dip, rake = finfo.rake, depth = finfo.depth;
	aziV.clear(); grtM.clear(); phtM.clear(); ampM.clear();
   //float azi[nazi], grT[nper][nazi], phT[nper][nazi], amp[nper][nazi];
	aziV = std::vector<float>( azi, azi+nazi );
	for( int iper=0; iper<perlst.size(); iper++ ) {
		float per = perlst[iper];
      grtM[per] = std::vector<float>( grT[iper], grT[iper]+nazi );
      phtM[per] = std::vector<float>( phT[iper], phT[iper]+nazi );
      ampM[per] = std::vector<float>( amp[iper], amp[iper]+nazi );
		
/*
      for(int iazi=0; iazi<nazi; iazi++) {
         if( azi[iazi] != iazi*2 )
	    throw ErrorRP::BadAzi(FuncName, std::to_string(azi[iazi])+" != 2*"+std::to_string(iazi));
         per_azi_pred[i][iazi] = AziData( azi[iazi], grT[i][iazi], phT[i][iazi], amp[i][iazi] );
      }
*/
   }
}


void RadPattern::GetPred( const float per, const float azi,
								  float& grt, float& pht, float& amp ) const {
	// check validities of period and azimuth
	auto Igrt = grtM.find(per);
	if( Igrt == grtM.end() )
		throw ErrorRP::BadParam(FuncName, "un-predicted period");
	int iazi = (int)(azi/dazi);
	if( iazi<0 || iazi >= nazi )
		throw ErrorRP::BadAzi( FuncName, "azi = "+std::to_string(azi) );
	// low and high azimuth
	float azil = iazi*dazi, azih = azil+dazi;
	float azifactor = (azi-azil) / dazi, ftmp1, ftmp2;
	// group delay
	ftmp1 = (Igrt->second)[iazi], ftmp2 = (Igrt->second)[iazi+1];
	grt = ftmp1 + (ftmp2 - ftmp1) * azifactor;
	// phase shift
	auto Ipht = phtM.find(per);
	ftmp1 = (Ipht->second)[iazi], ftmp2 = (Ipht->second)[iazi+1];
	pht = ftmp1 + (ftmp2 - ftmp1) * azifactor;
	// amp
	auto Iamp = ampM.find(per);
	ftmp1 = (Iamp->second)[iazi], ftmp2 = (Iamp->second)[iazi+1];
	amp = ftmp1 + (ftmp2 - ftmp1) * azifactor;

}
