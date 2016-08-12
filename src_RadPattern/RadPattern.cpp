#include "RadPattern.h"
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>

/* FORTRAN entrance */
const int nazi = RadPattern::nazi;
extern"C" {
   //void rad_pattern_r_(const float *strike, const float *dip, const float *rake, int *nper, const float *dper, 
   void rad_pattern_r_(const float mt[6], int *nper, const float *dper, 
							  const float *per, const float *eigH, const float *deigH, const float *eigV, const float *deigV, const float *camp, const float *wvn,
							  const float *perlst, int *nperlst, float *azi, float grT[][nazi], float phT[][nazi], float amp[][nazi]);
   //void rad_pattern_l_(const float *strike, const float *dip, const float *rake, int *nper, const float *dper,
   void rad_pattern_l_(const float mt[6], int *nper, const float *dper,
							  const float *per, const float *eigH, const float *deigH, const float *camp, const float *wvn,
							  const float *perlst, int *nperlst, float *azi, float grT[][nazi], float phT[][nazi], float amp[][nazi]);
}


/* con/destructors and operators */
RadPattern::RadPattern( const char type, const std::string& feigname )
   : type(type), er(feigname, 1, true) {
	if( type == 'N' ) return;
   if( type!='R' && type!='L' )
      throw ErrorRP::BadParam(FuncName, "unknown type = "+type);
	er.FillSD();
}

void RadPattern::SetModel( const char typein, const std::string& feigname ) {
   type = typein;	er.reLoad(feigname, 1, true); er.FillSD();
}

/* copy arrayin[nazi] into Vout with positions shifted by int(nazi/2) */
void RadPattern::ShiftCopy( std::vector<float>& Vout, const float* arrayin, const int nazi ) const {
   // check nazi
   if( nazi % 2 == 0 )
      throw ErrorRP::BadParam( FuncName, "unexpected nazi = " + std::to_string(nazi) );
   int nazio2 = nazi / 2;
   Vout.clear(); Vout.reserve(nazi);
   // shift by 180 degree
   Vout = std::vector<float>( arrayin+nazio2, arrayin+nazi-1 );
   Vout.insert( Vout.end(), arrayin, arrayin+nazio2+1 );
}

// computes moment tensor from strike, dip, and rake
std::array<float, 6> RadPattern::MomentTensor( float stk, float dip, float rak, const float M0 ) const {
	float deg2rad = M_PI/180.;
	stk *= deg2rad; dip *= deg2rad; rak *= deg2rad;
	float sins = sin(stk), coss = cos(stk), sin2s = sin(2.*stk), cos2s = cos(2.*stk);
	float sind = sin(dip), cosd = cos(dip), sin2d = sin(2.*dip), cos2d = cos(2.*dip);
	float sinr = sin(rak), cosr = cos(rak);
	std::array<float, 6> MT;
	MT[0] = -M0 * (sind*cosr*sin2s + sin2d*sinr*sins*sins);	// xx
	MT[1] =  M0 * (sind*cosr*sin2s - sin2d*sinr*coss*coss);	// yy
	MT[2] =  M0 * (sin2d*sinr);										// zz
	MT[3] =  M0 * (sind*cosr*cos2s + sin2d*sinr*sins*coss);	// xy
	MT[4] = -M0 * (cosd*cosr*coss + cos2d*sinr*sins);			// xz
	MT[5] = -M0 * (cosd*cosr*sins - cos2d*sinr*coss);			// yz
	return MT;
}

/* predict radpattern for rayleigh and love waves */
bool RadPattern::Predict( const std::array<ftype, 6>& MTi, const ftype depin, const ftype M0in, const std::vector<float>& perlst ) {

	// return if the requested new state is exactly the same as the one stored
	//if( stk==stkin && dip==dipin && rak==rakin &&
	if( MT[0]==MTi[0] && MT[1]==MTi[1] && MT[2]==MTi[2] && MT[3]==MTi[3] && MT[4]==MTi[4] && MT[5]==MTi[5] &&
		 dep==depin && M0==M0in && perlst.size()<=grtM.size() ) {
		bool allfound = true;
		for( const auto per : perlst )
			if( grtM.find(per) == grtM.end() ) {
				allfound = false;
				break;
			}
		if( allfound ) return false;	// not updated
	}

	// store current state;
	/*
	stk = stkin; dip = dipin;
	rak = rakin;
	*/
	MT = MTi; M0 = M0in; dep = depin;

	// extract source data at dep
	er.FillSDAtDep( dep );

	// TODO extract wvn and coef_amp
	campM.clear();
	for( int iper=0; iper<perlst.size(); iper++ ) {
		float per = perlst[iper];
		auto Iper = std::lower_bound(er.sd.per.begin(), er.sd.per.end(), per);
		if( Iper == er.sd.per.end() )
			throw ErrorRP::BadParam( FuncName, "non-exist period = " + std::to_string(per) );
		auto index = Iper - er.sd.per.begin();
		campM[per] = std::array<float, 2>{ er.sd.ac[index], er.sd.wvn[index] };
   }

   // run rad_pattern
   int nperlst = perlst.size();
	float azi[nazi], grT[nperlst][nazi], phT[nperlst][nazi], amp[nperlst][nazi];

   if( type == 'R' ) {
std::cerr<<"debugR: "<<MT[0]<<" "<<MT[1]<<" "<<MT[2]<<" "<<MT[3]<<" "<<MT[4]<<" "<<MT[5]<<std::endl;
      rad_pattern_r_( MT.data(), &(er.sd.nper), &(er.sd.dper), 
							 er.sd.per.data(), er.sd.eigH.data(), er.sd.deigH.data(), er.sd.eigV.data(), er.sd.deigV.data(), er.sd.ac.data(), er.sd.wvn.data(),
							 perlst.data(), &nperlst, azi, grT, phT, amp );
   } else if( type == 'L' ) {
std::cerr<<"debugL: "<<MT[0]<<" "<<MT[1]<<" "<<MT[2]<<" "<<MT[3]<<" "<<MT[4]<<" "<<MT[5]<<std::endl;
      rad_pattern_l_( MT.data(), &(er.sd.nper), &(er.sd.dper), 
							 er.sd.per.data(), er.sd.eigH.data(), er.sd.deigH.data(), er.sd.ac.data(), er.sd.wvn.data(),
							 perlst.data(), &nperlst, azi, grT, phT, amp );
		//for(int i=0; i<nazi; i++) std::cout<<azi[i]<<" "<<grT[0][i]<<" "<<phT[0][i]<<" "<<amp[0][i]<<std::endl;
   }

	// copy predictions into maps
	grtM.clear(); phtM.clear(); ampM.clear(); //campM.clear(); aziV.clear();
   //float azi[nazi], grT[nperlst][nazi], phT[nperlst][nazi], amp[nperlst][nazi];
	for( int iper=0; iper<perlst.size(); iper++ ) {
		float per = perlst[iper];
/*
		// copy as is
		auto &grV = grtM[per], &phV = phtM[per], &amV = ampM[per];
      grV = std::vector<float>( grT[iper], grT[iper]+nazi );
      phV = std::vector<float>( phT[iper], phT[iper]+nazi );
      amV = std::vector<float>( amp[iper], amp[iper]+nazi );
		// shift the phase by pi
		float T_pi = per * 0.5;
		auto pi_shifter = [&](float& val) {
			val += T_pi;
			if( val > T_pi ) val -= per;
		};
		std::for_each( phV.begin(), phV.end(), pi_shifter );
		// flip
		std::transform( grV.begin(), grV.end(), grV.begin(), std::negate<float>() );
		std::transform( phV.begin(), phV.end(), phV.begin(), std::negate<float>() );
*/
		// shift by 180 degree
		ShiftCopy( grtM[per], grT[iper], nazi );
		ShiftCopy( phtM[per], phT[iper], nazi );
		ShiftCopy( ampM[per], amp[iper], nazi ); for( auto &val : ampM[per] ) val *= M0;
		//campM[per] = std::array<float, 2>{ coef_amp[iper], wvn[iper] };
   }
	//aziV = std::vector<float>( azi, azi+nazi );

	// invalidate focal predictions with amplitudes < amp_avg * AmpValidPerc
	//const float NaN = Rimpl::NaN;
	for( int iper=0; iper<perlst.size(); iper++ ) {
		float per = perlst[iper];
		// determine min amplitude
		int Nvalid = 0; float Amin = 0.;
		const auto& ampV = ampM[per];
		for( const auto& amp : ampV ) Amin += amp;
		Amin *= (AmpValidPerc / ampV.size());
		// invalidate azimuths with small amplitudes
		auto& grtV = grtM[per];
		for(int iazi=0; iazi<nazi; iazi++)
			if( ampV[iazi] < Amin ) {
				int jazilow = iazi-InvalidateHwidth, jazihigh = iazi+InvalidateHwidth+1;
				// invalidate the range (jazilow - 0)
				int jazi = jazilow;
				if( jazi < 0 )
					for(; jazi<0; jazi++)
						grtV[jazi + nazi] = NaN;
				// invalidate the range (0 - 360)
				for(; jazi<jazihigh&&jazi<nazi; jazi++)
					grtV[jazi] = NaN;
				// invalidate the range (360 - jazihigh)
				if( jazihigh > nazi )
					for(; jazi<jazihigh; jazi++)
						grtV[jazi - nazi] = NaN;
			}
	}
	return true;	// updated!
}

std::array<float, 2> RadPattern::cAmp( const float per ) const {
	auto iter = campM.find(per);
	if( iter == campM.end() )
		throw ErrorRP::BadParam(FuncName, "un-predicted period");
	return (iter->second);
}

bool RadPattern::GetPred( const float per, const float azi,
								  float& grt, float& pht, float& amp,
								  const float Amul ) const {
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
	//const float NaN = Rimpl::NaN;
	// group delay
	ftmp1 = (Igrt->second)[iazi], ftmp2 = (Igrt->second)[iazi+1];
	if( ftmp1==NaN || ftmp2==NaN ) {
		grt = pht = amp = NaN;
		return false;
	}
	grt = ftmp1 + (ftmp2 - ftmp1) * azifactor;
	// phase shift
	auto Ipht = phtM.find(per);
	ftmp1 = (Ipht->second)[iazi], ftmp2 = (Ipht->second)[iazi+1];
	pht = ftmp1 + (ftmp2 - ftmp1) * azifactor;
	// amp
	auto Iamp = ampM.find(per);
	ftmp1 = (Iamp->second)[iazi], ftmp2 = (Iamp->second)[iazi+1];
	amp = Amul * (ftmp1 + (ftmp2 - ftmp1) * azifactor);

	return true;
}

bool RadPattern::GetPred( const float per, const float azi,
								  float& grt, float& pht, float& amp,
								  const float dis, const float alpha, const float recCAmp ) const {
	if( ! GetPred(per, azi, grt, pht, amp, 1.) ) return false;

	//if( U==NaN || J==NaN ) {	// compute amplitude at the source
	if( recCAmp == NaN ) {	// compute amplitude at the source
		// M0 = scalar seismic momentum
		// cAmp = source amp norm term: 1.0/( (phvel*grvel*I0) * sqrt(8 * pi) ) in sec^2/gram
		auto coefs = cAmp(per);	// cAmp and wavenumber
		//std::cerr<<"amp0 = "<<amp<<" ";
		amp *= 100. * coefs[0] / sqrt(coefs[1]);	// amplitude at 1km distance
		//std::cerr<<100.<<" "<<M0<<" "<<sqrt(coefs[0]*1.0e15/sqrt(8*M_PI))<<" "<<1.0e-6*sqrt(coefs[0]*1.0e15*sqrt(8*M_PI))/sqrt(coefs[1])*exp(-dis*alpha)<<" amp1 = "<<amp<<" ";
	} else { // propogate amp to the receiver if extra info is given
		// recCAmp = 1. / sqrt(J * U * C) in sqrt(sec^2/gram)
		// J = receiver mode energy integration (from eigen);
		// U = local group velocity at the receiver location
		// C = local phase velocity at the receiver location
		auto coefs = cAmp(per);	// cAmp and wavenumber
		amp *= 100. * sqrt( coefs[0] / (sqrt(8.*M_PI)) )
				 * recCAmp / sqrt(coefs[1]);		// receiver norm term * second part of propogation
	}
	if( dis>0 ) {	// attenuation and geometric spreading
		// dis = distance; alpha = average attenuation coeff
		amp /= sqrt(dis);
		if( alpha > 0 ) amp *= exp(-dis*alpha);		// anelastic propogation
		//std::cerr<<dis<<" "<<alpha<<" amp2 = "<<amp<<"  ";
	}

	return true;
}

void RadPattern::OutputPreds( const std::string& fname, const float norm_dis, const float Q ) const {
	if( grtM.size() == 0 ) return;

	std::ofstream fout( fname );
   if( ! fout ) throw ErrorRP::BadFile(FuncName, fname);

	for( auto &grtP : grtM ) {
		float per = grtP.first;
		float alpha = M_PI/(per*3.0*Q);
		auto &grtV = grtP.second;
		auto &phtV = phtM.at(per);
		auto &ampV = ampM.at(per);
		int iazi; float azi, grt, pht, amp;
		auto getPreds = norm_dis>0 ? std::function<bool()>([&](){ return GetPred(per, azi, grt, pht, amp, norm_dis, alpha); }) : 
							 [&](){ grt=grtV.at(iazi); if(grt==RadPattern::NaN)return false; pht=phtV.at(iazi); amp=ampV.at(iazi); return true; };
		for( iazi=0; iazi<nazi; iazi++ ) {
			azi = iazi*dazi;
			if( getPreds() ) fout<<azi<<" "<<grt<<" "<<pht<<" "<<amp<<" "<<per<<"\n";
		}
		fout<<"\n\n";
	}
}

RadPattern& RadPattern::operator-=(const RadPattern &rp2) {
	for( auto &pvPair : grtM ) {
		// check for consistency
		auto per = pvPair.first; auto size = pvPair.second.size();
		auto I2 = rp2.grtM.find(per);
		if( (I2->first)!=per || (I2->second).size()!=size ) {
			std::stringstream ss; ss<<"("<<per<<", "<<size<<") - ("<<(I2->first)<<", "<<(I2->second).size()<<")";
			throw ErrorRP::HeaderMismatch(FuncName, ss.str());
		}
		// locate group, phase, and amplitude vectors for current period
		auto grtA1 = pvPair.second.data(); auto grtA2 = I2->second.data();					// group
		auto phtA1 = phtM[per].data(); auto phtA2 = rp2.phtM.at(per).data();	// phase
		auto ampA1 = ampM[per].data(); auto ampA2 = rp2.ampM.at(per).data();	// amplitude
		for(int i=0; i<size; i++) {
			if( grtA1[i]==NaN ) continue;
			if( grtA2[i]==NaN ) { grtA1[i] = NaN; continue; }
			grtA1[i] -= grtA2[i]; phtA1[i] -= phtA2[i]; ampA1[i] -= ampA2[i];
		}
	}
	return *this;
}

RadPattern& RadPattern::operator*=(const float &mul) {
	for( auto &pvPair : grtM ) {
		// check for consistency
		auto grtA = pvPair.second.data();
		auto ampA = ampM.at(pvPair.first).data();
		for(int i=0; i<nazi; i++)
			if( grtA[i] != RadPattern::NaN ) ampA[i] *= mul;
	}
	M0 *= mul;
	return *this;
}

// re-scale (the amplitude of) *this to have the minimum amp chiSquare misfit to rp2
// sigmaV is a map (by period) of sigmaAs (fraction 0-1)
void RadPattern::NormBy( const RadPattern &rp2, const MA3& sigmasM ) {
	std::map<float, float> sigmaM;
	for( const auto &paPair : sigmasM ) sigmaM[paPair.first] = (paPair.second)[2];
	NormBy(rp2, sigmaM);
}
void RadPattern::NormBy( const RadPattern &rp2, const std::map<float,float>& sigmaM ) {
	if( grtM.empty() ) return;
	float a = 0., b = 0.;
	NormCoefs( rp2, sigmaM, a, b );
	//std::cerr<<"NormBy: "<<a<<" "<<b<<" "<<exp(b/a)<<std::endl;
	*this *= exp(b/a);
}
void NormAmps( RadPattern &rp1R, RadPattern &rp1L, const RadPattern &rp2R, const RadPattern &rp2L,
					const MA3 &sigmasMR, const MA3 &sigmasML ) {
	std::map<float, float> sigmaMR, sigmaML;
	for( const auto &paPair : sigmasMR ) sigmaMR[paPair.first] = (paPair.second)[2];
	for( const auto &paPair : sigmasML ) sigmaML[paPair.first] = (paPair.second)[2];
	float a = 0., b = 0.;
	rp1R.NormCoefs( rp2R, sigmaMR, a, b );
	rp1L.NormCoefs( rp2L, sigmaML, a, b );
	float mul = exp(b/a); rp1R *= mul; rp1L *= mul;
}

void RadPattern::NormCoefs( const RadPattern &rp2, const std::map<float,float>& sigmaM, float &a, float &b ) {
	for( auto &pvPair : grtM ) {
		// check for consistency
		auto per = pvPair.first; auto size = pvPair.second.size();
		auto I2 = rp2.grtM.find(per);
		if( I2==rp2.grtM.end() || (I2->second).size()!=size ) {
			std::stringstream ss; ss<<"("<<per<<", "<<size<<") - ("<<(I2->first)<<", "<<(I2->second).size()<<")";
			throw ErrorRP::HeaderMismatch(FuncName, ss.str());
		}
		// accumulate a and b
		float weightA = -log(1.-sigmaM.at(per)); weightA = 1. / (weightA*weightA); 
		// locate group and amplitude vectors for current period
		auto grtA1 = pvPair.second.data(); auto grtA2 = I2->second.data();	// group
		auto ampA1 = ampM[per].data(); auto ampA2 = rp2.ampM.at(per).data();	// amplitude
		int n = 0; float l1 = 0., l2 = 0.;
		for(int i=0; i<size; ++i) {
			if( grtA1[i]==RadPattern::NaN || grtA2[i]==RadPattern::NaN ) continue;
			l1 += log(ampA1[i]); l2 += log(ampA2[i]);	++n;
		}
		a += n * weightA;	b += (l2-l1)*weightA; 
	}
}

// sigmasV is a map (by period) of sigma triples: sigmaG (sec), sigmaP (sec), sigmaA (fraction 0-1)
inline int nint(const float& val) { return (int)floor(val+0.5); }
std::array<float,4> RadPattern::chiSquare(const RadPattern &rp2, const MA3 &sigmasM, const bool normalizeA) {
	if( grtM.empty() ) return {NaN, NaN, NaN, NaN};
	// normalize the amplitudes of rp2in if requested
	//RadPattern rp2C;
	if( normalizeA ) { NormBy(rp2, sigmasM); }
	//auto &rp2 = normalizeA ? rp2C : rp2in;

	int n = 0; float chiSG = 0., chiSP = 0., chiSA = 0.;
	for( auto &pvPair : grtM ) {
		// check for consistency
		auto per = pvPair.first; auto size = pvPair.second.size();
		auto I2 = rp2.grtM.find(per);
		if( (I2->first)!=per || (I2->second).size()!=size ) {
			std::stringstream ss; ss<<"("<<per<<", "<<size<<") - ("<<(I2->first)<<", "<<(I2->second).size()<<")";
			throw ErrorRP::HeaderMismatch(FuncName, ss.str());
		}
		// locate group, phase, and amplitude vectors for current period
		auto grtA1 = pvPair.second.data(); auto grtA2 = I2->second.data();	// group
		auto phtA1 = phtM.at(per).data(); auto phtA2 = rp2.phtM.at(per).data();	// phase
		auto ampA1 = ampM.at(per).data(); auto ampA2 = rp2.ampM.at(per).data();	// amplitude
		// compute rms
		float rmsG = 0., rmsP = 0., rmsA = 0.;
		for(int i=0; i<size; i++) {
			if( grtA1[i]==RadPattern::NaN || grtA2[i]==RadPattern::NaN ) continue;
			float gdiff = grtA1[i] - grtA2[i]; rmsG += gdiff*gdiff;
			float pdiff = phtA1[i] - phtA2[i]; pdiff -= nint(pdiff/per)*per; rmsP += pdiff*pdiff;
			float Adiff = log(ampA1[i]/ampA2[i]); rmsA += Adiff*Adiff;
			++n;
		}
		// chi square formulas
		auto &sigmasA = sigmasM.at(per);	// get sigmas for the current period
		float sigmaA = -log(1.-sigmasA[2]);
		chiSG += rmsG / (sigmasA[0]*sigmasA[0]);
		chiSP += rmsP / (sigmasA[1]*sigmasA[1]);
		chiSA += rmsA / (sigmaA*sigmaA);
	}

	return std::array<float, 4>{chiSG, chiSP, chiSA, (float)n};
}

std::array<float,4> chiSquare( RadPattern &rp1R, RadPattern &rp1L, const RadPattern &rp2R, const RadPattern &rp2L,
										  const MA3 &sigmasMR, const MA3 &sigmasML, const bool normalizeA ) {
	if(normalizeA) NormAmps( rp1R, rp1L, rp2R, rp2L, sigmasMR, sigmasML );
	auto resAR = rp1R.chiSquare( rp2R, sigmasMR, false );
	auto resAL = rp1L.chiSquare( rp2L, sigmasML, false );
	int n = 0; float chiSG = 0., chiSP = 0., chiSA = 0.;
	if( resAR[0] != RadPattern::NaN ) { chiSG += resAR[0]; chiSP += resAR[1]; chiSA += resAR[2]; n += resAR[3]; }
	if( resAL[0] != RadPattern::NaN ) { chiSG += resAL[0]; chiSP += resAL[1]; chiSA += resAL[2]; n += resAL[3]; }
	return std::array<float, 4>{chiSG, chiSP, chiSA, (float)n};
}

