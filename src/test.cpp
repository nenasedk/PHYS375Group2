#include "AdaptSolve.hh"
#include <vector>

using namespace std;

void DYDX(double x, vector<double>& y, vector<double>& dydx){
  dydx.at(0) = x;
}

int main(){
  vector<double> y = vector<double>(2,0.0);
  double start = 1e-10;
  double end = 10.0;
  AdaptSolve a = AdaptSolve();
  a.RKSolve(y,2,start,end,DYDX);
  cout << a.kount << endl;
  for(int i=0;i<a.kount;i++){
    cout << "x: " << a.xp.at(i) << " y: " << a.yp.at(0).at(i) << endl;
  }

  return 0;
}
  
