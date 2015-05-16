#include <iostream>
#include <fstream>

int main() {
	std::ifstream fin("test.txt");
	std::string line;
	std::getline(fin, line);
	float lon, lat, data;
	std::cout<< sscanf(line.c_str(), "%f %f %f", &lon, &lat, &data) <<"\n";
	std::cout<<lon<<" "<<lat<<" "<<data<<"\n";
}
