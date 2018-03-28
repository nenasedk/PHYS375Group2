#include "Star.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include "math.h"
using namespace std;

int main(){
	
for(int loop = 1; loop < 101; loop++){
	double unitlessTemp = 1;
	double unitlessDens = loop;
	double X = 0.734;
	double Y = 0.250;
	double Z = 0.016;
	double mu = pow((2.0*X + 0.75*Y + 0.5*Z),-1);
	Star(unitlessDens, unitlessTemp, X, Y, Z, mu);
	std::string filename("MSStar");
	std::stringstream fileNameNum;
	fileNameNum << filename << loop << ".txt";
	std::string fileNameResult = fileNameNum.str();
	
	//cout << fileNameResult << endl;
	ofstream myfile (fileNameResult);
	if (myfile.isopen()){
		myfile << _Rad << _Dens << _Temp << _Mass << _Lum << _OptD << endl;
	}
	
	myfile.close();
	
	}
	return 0;
}
