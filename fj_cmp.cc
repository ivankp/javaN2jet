#include <iostream>
#include <iomanip>
#include <sstream>
#include <array>
#include <vector>
#include <functional>
#include <cstring>
#include <cstdio>
#include <sys/time.h>

#include <fastjet/ClusterSequence.hh>

using namespace std;
using namespace fastjet;

#include <TRandom3.h>

using namespace std;

string exec(const char* cmd) {
  FILE* pipe = popen(cmd, "r");
  if (!pipe) return "ERROR";
  char buffer[128];
  std::string result = "";
  while(!feof(pipe)) {
  	if(fgets(buffer, 128, pipe) != NULL)
  		result += buffer;
  }
  pclose(pipe);
  return result;
}

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

  stringstream cmd_n2;
  cmd_n2 << "java test2 " << argv[1] << ' ' << argv[2];

  while (true) {
    // Generate 2 random 4-vectors
    timeval tv;
    gettimeofday(&tv, nullptr);
    TRandom3 r(tv.tv_usec);
    stringstream ps;
    ps << fixed << scientific << setprecision(8);
    for (int i=0; i<8; ++i) ps << ' ' << r.Rndm();
    const string ps_str( ps.str() );

    // FastJet **********************************************
    const double R = atof(argv[2]);
    JetAlgorithm jalg;
    function<double(const PseudoJet& a, const PseudoJet& b)> dij;
    function<double(const PseudoJet& a)> diB;

    if (!strcmp(argv[1],"antikt")) {
      jalg = antikt_algorithm;
      dij = [R](const PseudoJet& a, const PseudoJet& b){
        return min(1./a.kt2(),1./b.kt2())*a.delta_R(b)/R;
      };
      diB = [](const PseudoJet& a){ return 1./a.kt2(); };
    } else if (!strcmp(argv[1],"kt")) {
      jalg = kt_algorithm;
      dij = [R](const PseudoJet& a, const PseudoJet& b){
        return min(a.kt2(),b.kt2())*a.delta_R(b)/R;
      };
      diB = [](const PseudoJet& a){ return a.kt2(); };
    } else if (!strcmp(argv[1],"cambridge")) {
      jalg = cambridge_algorithm;
      dij = [R](const PseudoJet& a, const PseudoJet& b){
        return a.delta_R(b)/R;
      };
      diB = [](const PseudoJet& a){ return 1.; };
    } else {
      cerr << "Undefined algorithm " << argv[1] << endl;
      return 1;
    }

    // collect particles
    vector<PseudoJet> particles;
    while (ps >> particles) { };

    // define ClusterSequence
    ClusterSequence seq(particles, JetDefinition(jalg,R), false);

    // cluster jets
    vector<PseudoJet> jets = seq.inclusive_jets();

    // sort jets by pT
    sort( jets.begin(), jets.end(),
      [](const PseudoJet& i, const PseudoJet& j)
        { return i.pt() > j.pt(); }
    );

    const string out_fj = [&jets](){
      stringstream ss_fj;
      ss_fj << fixed << scientific << setprecision(8);
      for (auto& jet : jets) ss_fj << jet.pt() << ' ';
      return ss_fj.str();
    }();

    // n2jet_java *******************************************

    const string out_n2 = exec((cmd_n2.str()+ps_str).c_str());

    // Print output *****************************************

    if (out_n2!=out_fj) {
      cout << "FJ: " << out_fj << '\n'
           << "N2: " << out_n2 << '\n'
           << "p1: " << ps_str.substr(1,59)
           << " diB = " << diB(particles[0]) << '\n'
           << "p2: " << ps_str.substr(61)
           << " diB = " << diB(particles[1]) << '\n'
           << "dij = " << dij(particles[0],particles[1]) << '\n'
           << endl;
    }
  }

  return 0;
}
