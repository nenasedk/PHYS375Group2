#include "AdaptSolve.hh"
#include <vector>

using namespace std;

void DYDX(double x, vector<double>& y, vector<double>& dydx){
  dydx.at(0) = x;
}
void DYDX2(double x, vector<double>& y, vector<double>& dydx){
  dydx.at(0) = 1;
}

int main(){
  vector<double> y = vector<double>(3,0.0);
  void (* ds[])(double, vector<double>&, vector<double>&) = {DYDX,DYDX2};
  double start = 1e-10;
  double end = 10.0;
  AdaptSolve a = AdaptSolve();
  a.RKSolve(y,3,start,end,&ds);
  cout << a.kount << endl;
  for(int i=0;i<a.kount;i++){
    cout << "x: " << a.xp.at(i) << " y: " << a.yp.at(0).at(i) << endl;
  }

  return 0;
}
  
