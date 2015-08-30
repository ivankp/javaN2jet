#include <iostream>
#include <vector>
#include <cmath>

#include <fastjet/ClusterSequence.hh>

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
using namespace fastjet;

// read 4-momentum from input stream
istream& operator>>(istream &in, vector<PseudoJet>& p) {
  double px, py, pz, E;
  in >> px >> py >> pz >> E;
  if (in) p.emplace_back(px,py,pz,E);
  return in;
}
ostream& operator<<(ostream &out, const PseudoJet& p) {
  return out << p.px() << ' ' << p.py() << ' '
             << p.pz() << ' ' << p.E ();
}

int main(int argc, char **argv)
{
  if (argc!=3) {
    cout << "Usage: " << argv[0] << " algorithm radius" << endl;
    return 1;
  }

  // collect particles
  vector<PseudoJet> p;
  while (cin >> p) { };

  test(p[0])
  test(p[1])

  test(p[0].rap())
  test((0.5*log((p[0].E()+p[0].pz())/(p[0].E()-p[0].pz()))))
  test(p[1].rap())
  test((0.5*log((p[1].E()+p[1].pz())/(p[1].E()-p[1].pz()))))
  test(p[0].delta_R(p[1]))

  return 0;
}
