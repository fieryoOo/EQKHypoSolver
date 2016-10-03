#include "RadPattern.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>


inline int nint(float datain) { return (int)floor(datain+0.5); }

int main( int argc, char* argv[] ) {
   if( argc!=10 && argc!=13 ) {
      std::cerr<<"Usage: "<<argv[0]<<" [R/L] [eigen_file (.R/L)] [per_lst] [depth] [M0]\n"
					<<"option 1: [strike] [dip] [rake] [out_name]\n"
					<<"option 2: [Mxx] [Myy] [Mzz] [Mxy] [Mxz] [Myz] [out_name]"<<std::endl;
					//<<"[out_name] [norm_dis (optional)] [Q (optional)]"<<std::endl;
      exit(-1);
   }

   /* read in type */
   char type = argv[1][0];
   if( type != 'R' && type != 'L' ) {
      std::cerr<<"Unknown type: "<<type<<std::endl;
      exit(0);
   }

   /* read in per.lst */
   std::vector<float> perlst;
   std::ifstream fin(argv[3]);
   if( ! fin ) {
      std::cerr<<"Error(main): Cannot read from file "<<argv[3]<<std::endl;
      exit(0);
   }
   for(std::string line; std::getline(fin, line); ) {
      float pertmp;
      sscanf(line.c_str(), "%f", &pertmp);
      perlst.push_back(pertmp);
   }
   fin.close();
   std::cout<<"### "<<perlst.size()<<" periods read in. ###"<<std::endl;

   RadPattern rp( type, argv[2] );

	const float dep = atof(argv[4]), M0 = atof(argv[5]);
	std::string outname( argc==10 ? argv[9] : argv[12] );
	std::cout<<"### Input Focal info = ("<<dep<<" "<<M0<<" ";
	if( argc == 10 ) {
		const float stk = atof(argv[6]), dip = atof(argv[7]), rak = atof(argv[8]);
	   std::cout<<stk<<" "<<dip<<" "<<rak<<"). ###"<<std::endl;
		rp.Predict( stk, dip, rak, dep, M0, perlst );
	} else {
		float m1 = atof(argv[6]), m2 = atof(argv[7]), m3 = atof(argv[8]),m4 = atof(argv[9]),m5 = atof(argv[10]),m6 = atof(argv[11]);
	   std::cout<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<" "<<m5<<" "<<m6<<"). ###"<<std::endl;
		rp.Predict( std::array<float, 6>{m1,m2,m3,m4,m5,m6}, dep, M0, perlst );
	}

	float norm_dis = 1, Q = 100;
	//if( argc == 12 ) { norm_dis = atof(argv[10]); Q = atof(argv[11]); }
	//rp.CorrectPhaseC();
	rp.CorrectPhase();
	rp.OutputPreds( outname, norm_dis, Q );

   return 0;
}
