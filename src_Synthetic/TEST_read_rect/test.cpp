#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <stdexcept>

#ifdef __cplusplus
extern"C" {
#endif
      void read_rect_model_(char fname[255], int* nmod, double* p, int* ierr);

      //void read_rect_model_new_(char* fmod_buff, int* buffsize, int* nmod, double* p, int* ierr);
      void read_rect_model_new_(char fname[255], int* nmod, double* p, int* ierr);

#ifdef __cplusplus
}
#endif


class fstring : public std::string {
public:
	fstring() : std::string() {}
	fstring( const char* strin ) : std::string(strin) {}
	fstring( const std::string& strin ) : std::string(strin) {}
	fstring( const fstring& fstr2 ) : std::string(fstr2), pfstring(nullptr) {}
	~fstring() { clearfstring(); }

	char* f_str() const { return f_str( size()+1 ); }
	char* f_str( const size_t fsize ) const {
		if( fsize < size()+1 )
			throw std::runtime_error("Error(fstring::f_str): input size too small");
		clearfstring();
		pfstring = new char[fsize];
		sprintf(pfstring, "%s", c_str());
		std::fill(pfstring+size(), pfstring+fsize, ' ');
		return pfstring;
	}

	void assignf( const char* strin, const size_t strsize ) {
		int i=strsize-1;
		for(; i>=0; i--)
			if( strin[i] != ' ' ) break;
		this->assign( strin, i+1 );
	}

private:
	mutable char* pfstring = nullptr;
	void clearfstring() const { if(pfstring!=nullptr) delete [] pfstring; pfstring=nullptr; }
};



int main() {
	fstring name_fmodel("/work1/tianye/Syndat-1.1/data/WUSmap.25.bin");
	int ierr, nmod; double per;
	per = -123.; ierr = -123;
	nmod = 0; read_rect_model_(name_fmodel.f_str(255), &nmod, &per, &ierr);
	std::cerr<<" In main 1: "<<per<<" "<<ierr<<std::endl;
	nmod = 1; read_rect_model_(name_fmodel.f_str(255), &nmod, &per, &ierr);
	std::cerr<<" In main 2: "<<per<<" "<<ierr<<std::endl;
	nmod = -1; read_rect_model_(name_fmodel.f_str(255), &nmod, &per, &ierr);
	std::cerr<<" In main 3: "<<per<<" "<<ierr<<std::endl;

	per = -123.; ierr = -123;
	nmod = 0; read_rect_model_new_(name_fmodel.f_str(255), &nmod, &per, &ierr);
	std::cerr<<" In main 4: "<<per<<" "<<ierr<<std::endl;
	nmod = 1; read_rect_model_new_(name_fmodel.f_str(255), &nmod, &per, &ierr);
	std::cerr<<" In main 5: "<<per<<" "<<ierr<<std::endl;
	nmod = -1; read_rect_model_new_(name_fmodel.f_str(255), &nmod, &per, &ierr);
	std::cerr<<" In main 6: "<<per<<" "<<ierr<<std::endl;
	
	return 0;
}
