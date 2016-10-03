#include <iostream>
#include <fstream>
#include <complex>
#include <array>

class Matrix {
public:
	Matrix(double mxx, double myy, double mzz, double mxy=0, double mxz=0, double myz=0)
		: mxx(mxx), mxy(mxy), mxz(mxz), myx(mxy), myy(myy), myz(myz), mzx(mxz), mzy(myz), mzz(mzz) {}
	Matrix(double mxx, double mxy, double mxz, double myx, double myy, double myz, double mzx, double mzy, double mzz)
		: mxx(mxx), mxy(mxy), mxz(mxz), myx(myx), myy(myy), myz(myz), mzx(mzx), mzy(mzy), mzz(mzz) {}
	Matrix(std::istream &sin) { sin >> mxx >> mxy >> mxz >> myx >> myy >> myz >> mzx >> mzy >> mzz; }

	std::array<std::complex<double>, 3> EigenValues() {
		// det coefs
		double a2 = -(mxx + myy + mzz);
		double a1 = -(mxy*myx + myz*mzy + mxz*mzx - mxx*myy -mxx*mzz - myy*mzz);
		double a0 = -(mxx*myy*mzz + mxy*myz*mzx + mxz*myx*mzy - mxx*myz*mzy -  myy*mxz*mzx - mzz*mxy*myx);
		std::cout<<"Solving x^3 + "<<a2<<"*x^2 + "<<a1<<"*x + "<<a0<<" = 0"<<std::endl;

		// cubic formula
		double Q = (3*a1-a2*a2) / 9;
		double R = (9*a2*a1-27*a0-2*a2*a2*a2) / 54;
		std::complex<double> sqrtD = std::sqrt(std::complex<double>(Q*Q*Q + R*R));
		std::complex<double> S = std::pow( R + sqrtD, 1./3. );
		std::complex<double> T = std::pow( R - sqrtD, 1./3. );

		using namespace std::literals::complex_literals;
		auto x1 = -a2/3. + (S+T);
		auto x2 = -a2/3. - 0.5*(S+T) + 0.5i*sqrt(3.)*(S-T);
		auto x3 = -a2/3. - 0.5*(S+T) - 0.5i*sqrt(3.)*(S-T);

		return {x1,x2,x3};
	}

	Matrix &operator+=(const Matrix &m2) {
		mxx += m2.mxx; mxy += m2.mxy;	mxz += m2.mxz; 
		myx += m2.myx; myy += m2.myy; myz += m2.myz;
		mzx += m2.mzx; mzy += m2.mzy;	mzz += m2.mzz;
	}

	friend std::ostream &operator<<(std::ostream &o, const Matrix &m) {
		o << m.mxx <<" "<< m.mxy <<" "<< m.mxz <<"\n"<< m.myx <<" "<< m.myy <<" "<< m.myz <<"\n"<< m.mzx <<" "<< m.mzy <<" "<< m.mzz;
		return o;
	}

private:
	double mxx, mxy, mxz, myx, myy, myz, mzx, mzy, mzz;
};

int main(int argc, char *argv[]) {
	if( argc != 2 ) {
		std::cerr<<"Usage: "<<argv[0]<<" [fmatrix(3x3)]"<<std::endl;
		return -1;
	}

	// input
	std::ifstream fin(argv[1]);
	if( ! fin ) {
		std::cerr<<"IO failed on file "<<argv[1]<<std::endl;
		return -2;
	}
	Matrix m(fin); std::cout<<"Input matrix:\n"<<m<<std::endl;

	auto ev = m.EigenValues();
	std::cerr<<ev[0]<<" "<<ev[1]<<" "<<ev[2]<<std::endl;
	
	return 0;
}
