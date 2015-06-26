#include <iostream>
#include <string>
#include <sstream>

int main() {
	std::string path("/home/tianye/EQKHypoSolver/README");
	std::stringstream ss(path);
	std::string token;
	while( std::getline(ss, token, '/') )
		std::cerr<<token<<std::endl;
	std::cerr<<"final token: "<<token<<"\n";
}
