#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>

#include <fastjet/ClusterSequence.hh>

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using namespace std;
using namespace fastjet;

// read 4-momentum from input stream
istream& operator>>(istream &in, vector<PseudoJet>& p) {
  double px, py, pz, E;
  in >> px >> py >> pz >> E;
  if (in) {
    p.emplace_back(px,py,pz,E);
    p.back().set_user_index(p.size()-1);
  }
  return in;
}
ostream& operator<<(ostream &out, const PseudoJet& p) {
  return out << p.user_index() << ": "
             << p.rap() << ' ' << p.phi() << ' ' << p.pt();
}

int main(int argc, char **argv)
{
  if (argc!=3) {
    cout << "Usage: " << argv[0] << " algorithm radius" << endl;
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

  // collect particles
  vector<PseudoJet> particles;
  while (cin >> particles) { };
  
  for (auto& p : particles) cout << p << endl;

  // define ClusterSequence
  ClusterSequence seq(particles, jdef, false);

  // cluster jets
  vector<PseudoJet> jets = seq.inclusive_jets();

  // sort jets by pT
  sort( jets.begin(), jets.end(),
    [](const PseudoJet& i, const PseudoJet& j){ return i.pt() > j.pt(); }
  );

  for (auto& jet : jets) {
    cout << jet << endl;
    for (auto& c : jet.constituents())
      cout <<'\t'<< c << endl;
  }

  return 0;
}
