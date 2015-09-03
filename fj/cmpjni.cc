#include <iostream>
#include <iomanip>
#include <sstream>
#include <array>
#include <vector>
#include <functional>
#include <chrono>
#include <random>
#include <future>
#include <cstring>
#include <cstdio>
#include <cmath>

#include <fastjet/ClusterSequence.hh>

#include <jni.h>

using namespace std;
using namespace fastjet;

int main(int argc, char **argv)
{
  if (argc!=3 && argc!=4) {
    cout << "Usage: " << argv[0] << " algorithm radius N" << endl;
    return 1;
  }

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

  JetDefinition jdef(jalg,R);

  // Setup JNI
  JavaVM *jvm;  // Pointer to the JVM (Java Virtual Machine)
  JNIEnv *jenv; // Pointer to native interface
  //================== prepare loading of Java VM ============================
  JavaVMInitArgs vm_args;                        // Initialization arguments
  JavaVMOption* options = new JavaVMOption[1];   // JVM invocation options
  options[0].optionString = (char*)"-Djava.class.path=.."; // dir of .class
  vm_args.version = JNI_VERSION_1_6;             // minimum Java version
  vm_args.nOptions = 1;                          // number of options
  vm_args.options = options;
  vm_args.ignoreUnrecognized = false; // invalid options make the JVM init fail
  //=============== load and initialize Java VM and JNI interface =============
  if (JNI_CreateJavaVM(&jvm, (void**)&jenv, &vm_args) != JNI_OK) {
    cerr << "JVM not loaded" << endl;
    return 1;
  }
  delete options; // no longer need the initialisation options.
  //=============== Display JVM version =======================================
  cout << "JVM load succeeded: Version ";
  jint ver = jenv->GetVersion();
  cout << ((ver>>16)&0x0f) << "."<<(ver&0x0f) << endl;

  jclass jcl = jenv->FindClass("Test");
  if (!jcl) {
    cerr << "Java class not found" << endl;
    return 1;
  }
  jmethodID jfcn = jenv->GetStaticMethodID(
    jcl, "cluster", "(Ljava/lang/String;D[D)[D"
  );
  if (!jfcn) {
    cerr << "Java method not found" << endl;
    return 1;
  }

  // mersenne twister random number generator
  mt19937 gen(chrono::system_clock::now().time_since_epoch().count());
  uniform_real_distribution<double> dist(0.0,1.0);

  long long event = 0;
  size_t np = (argc>3 ? atoi(argv[3]) : 0);
  double px, py, pz, E, *_p = (np ? new double[np*4] : nullptr);
  vector<PseudoJet> pp;
  while (true) {
    size_t ip = 0;

    // Generate random physical 4-vectors
    if (np) {
      pp.clear();
      pp.reserve(np);
      for (size_t i=0, n=np; i<n; ++i) {
        double m, pt,
               eta = 10.*(acos(1.-2.*dist(gen))/M_PI-0.5),
               phi = 2.*M_PI*dist(gen);
        while (!isfinite(pt = 10.-150.*log(1.-dist(gen)))) { }
        while (!isfinite(m  =     -20.*log(1.-dist(gen)))) { }

        _p[ip++] = px = pt*cos(phi);
        _p[ip++] = py = pt*sin(phi);
        _p[ip++] = pz = pt*sinh(eta);
        _p[ip++] = E  = sqrt(px*px+py*py+pz*pz+m*m);

        pp.emplace_back(px,py,pz,E);
      }
    } else {
      while (cin >> px >> py >> pz >> E) pp.emplace_back(px,py,pz,E);
      for (const PseudoJet& p : pp) {
        _p[ip++] = p.px();
        _p[ip++] = p.py();
        _p[ip++] = p.E ();
        _p[ip++] = p.pz();
      }
      np = pp.size();
    }

    // FastJet **********************************************

    auto fut_fj = async(launch::async, [&jdef,&pp](){
      // cluster jets
      vector<PseudoJet> jets =
        ClusterSequence(pp, jdef, false).inclusive_jets();

      // sort jets by pT
      sort( jets.begin(), jets.end(),
        [](const PseudoJet& i, const PseudoJet& j){ return i.pt() > j.pt(); }
      );

      // return jets' pT
      size_t njets = jets.size();
      vector<double> pts(njets);
      for (size_t i=0; i<njets; ++i) pts[i] = jets[i].pt();
      return pts;
    });

    // n2jet_java *******************************************

    // auto fut_n2 = async(launch::async, [=]() {
      jdoubleArray jpp = jenv->NewDoubleArray(np*4);
      jenv->SetDoubleArrayRegion(jpp,0,np*4,_p);
      jdoubleArray jpts = (jdoubleArray) jenv->CallStaticObjectMethod(
        jcl, jfcn, jenv->NewStringUTF(argv[1]), (jdouble)R, jpp);

      jsize size = jenv->GetArrayLength(jpts);
      vector<double> out_n2(size);
      jenv->GetDoubleArrayRegion(jpts,0,size,out_n2.data());

    //   return pts;
    // });

    // Print output *****************************************

    const vector<double> out_fj = fut_fj.get();
    // const vector<double> out_n2 = fut_n2.get();

    if (out_n2!=out_fj || argc==3) {
      const auto flags = cout.flags();
      cout << right << fixed << scientific << setprecision(8);
      cout << "\nFJ:";
      for (double pt : out_fj) cout <<' '<< pt;
      cout << endl;
      cout << "N2:";
      for (double pt : out_n2) cout <<' '<< pt;
      cout << endl;
      cout.flags(flags);

      for (size_t i=0; i<np; ++i) {
        cout << "p" << i << ": ";
        cout << right << fixed << scientific << setprecision(8);
        cout << pp[i].px() << ' ' << pp[i].py() << ' '
             << pp[i].pz() << ' ' << pp[i].E ();
        cout.flags(flags);
        cout << " diB = " << diB(pp[i]) << endl;
      }
      if (np<21)
        for (size_t i=1; i<np; ++i)
          for (size_t j=0; j<i; ++j)
            cout << "d" << setw(3) << i << setw(3) << j
                 << " = " << dij(pp[i],pp[j]) << endl;
      cout << endl;
    } else {
      for (int i=0;i<9;++i) cout << '\b';
      cout << setw(9) << event;
      cout.flush();
    }

    ++event;

    if (argc==3) break;
  }

  jvm->DestroyJavaVM();

  return 0;
}
