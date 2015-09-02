#include <iostream>
#include <iomanip>
#include <chrono>
#include <random>
#include <cstring>
#include <cmath>

#include <fastjet/ClusterSequence.hh>

using namespace std;
using namespace fastjet;

class running_stat {
  int n;
  double oldM, newM, oldS, newS;

public:
  running_stat(): n(0), oldM(0.), newM(0.), oldS(0.), newS(0.) { }
  void clear() noexcept { n=0; }
  void push(double x) noexcept {
    n++;
    // See Knuth TAOCP vol 2, 3rd edition, page 232
    if (n == 1) {
      oldM = newM = x;
      oldS = 0.;
    } else {
      newM = oldM + (x - oldM)/n;
      newS = oldS + (x - oldM)*(x - newM);
      // set up for next iteration
      oldM = newM;
      oldS = newS;
    }
  }
  int      num() const noexcept { return n; }
  double  mean() const noexcept { return (n > 0) ? newM : 0.; }
  double   var() const noexcept { return ( (n > 1) ? newS/(n - 1) : 0. ); }
  double stdev() const noexcept { return sqrt( var() ); }
};

int main(int argc, char **argv) {
  // check arguments
  if (argc!=5) {
    cout << "Usage: "<<argv[0]<<" [kt,antikt,cambridge] R Np Nev" << endl;
    return 1;
  }
  
  JetDefinition jdef(
    ( strcmp(argv[1],"antikt")
    ? ( strcmp(argv[1],"kt")
      ? ( strcmp(argv[1],"cambridge")
        ? throw runtime_error(string("Undefined algorithm ")+argv[1])
        : cambridge_algorithm )
      : kt_algorithm )
    : antikt_algorithm ), atof(argv[2])
  );

  // mersenne twister random number generator
  mt19937 gen(chrono::system_clock::now().time_since_epoch().count());
  uniform_real_distribution<double> dist(0.0,1.0);

  int np  = atoi(argv[3]);
  int nev = atoi(argv[4]);

  running_stat stat;

  for (int ev=0; ev<nev; ++ev) {
    vector<PseudoJet> pp;
    pp.reserve(np);

    for (int p=0; p<np; ++p) {
      double m, pt,
             eta = 10.*(acos(1.-2.*dist(gen))/M_PI-0.5),
             phi = 2.*M_PI*dist(gen);
      while (!isfinite(pt = 10.-150.*log(1.-dist(gen)))) { }
      while (!isfinite(m  =     -20.*log(1.-dist(gen)))) { }

      double px = pt*cos(phi),
             py = pt*sin(phi),
             pz = pt*sinh(eta),
             E  = sqrt(px*px+py*py+pz*pz+m*m);

      pp.emplace_back(px,py,pz,E);
    }

    using namespace std::chrono;
    auto t1 = high_resolution_clock::now();
    vector<PseudoJet> jets = ClusterSequence(pp,jdef,false).inclusive_jets();
    stat.push(duration_cast<nanoseconds>(high_resolution_clock::now() - t1).count()/1000.);
  }

  cout << stat.mean() << " ± " << stat.stdev() << " μs" << endl;
}
