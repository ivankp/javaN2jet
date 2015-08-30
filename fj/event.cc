#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <sys/time.h>

#include <fastjet/ClusterSequence.hh>

using namespace std;
using namespace fastjet;

struct Timer {
  timeval timer[2];

  timeval start() {
    gettimeofday(&this->timer[0], NULL);
    return this->timer[0];
  }
  timeval stop() {
    gettimeofday(&this->timer[1], NULL);
    return this->timer[1];
  }
  int duration() const {
    int secs(this->timer[1].tv_sec - this->timer[0].tv_sec);
    int usecs(this->timer[1].tv_usec - this->timer[0].tv_usec);
    if(usecs < 0) {
      --secs;
      usecs += 1000000;
    }
    return static_cast<int>(secs * 1000 + usecs / 1000.0 + 0.5);
  }
};

// read 4-momentum from input stream
istream& operator>>(istream &in, vector<PseudoJet>& p) {
  double px, py, pz, E;
  in >> px >> py >> pz >> E;
  if (in) p.emplace_back(px,py,pz,E);
  return in;
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

  Timer t;
  t.start();

  // define ClusterSequence
  ClusterSequence seq(particles, jdef, false);

  // cluster jets
  vector<PseudoJet> jets = seq.inclusive_jets();

  // sort jets by pT
  sort( jets.begin(), jets.end(),
    [](const PseudoJet& i, const PseudoJet& j){ return i.pt() > j.pt(); }
  );

  t.stop();
  cout << "FJ run time: " << t.duration() << " ms" << endl;
  cout << fixed << scientific << setprecision(8);
  for (auto& jet : jets) cout << jet.pt() << endl;

  return 0;
}
