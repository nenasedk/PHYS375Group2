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
		for (int i = 0; i < _Rad.size(); i++){
			myfile << _Rad(i) << "," << _Dens(i) << "," << _Temp(i) << "," << _Mass(i) << "," << _Lum(i) << "," << _OptD(i) << endl;	
		}
	}
	
	myfile.close();
	
	}
	return 0;
}
