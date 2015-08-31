#include <iostream>
#include <iomanip>

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
  // collect particles
  vector<PseudoJet> particles;
  while (cin >> particles) { };
  
  PseudoJet sum(0,0,0,0);
  for (auto& p : particles) sum += p;

  cout << fixed << scientific << setprecision(8) << sum << endl;

  return 0;
}
