import java.util.*;
import java.io.*;

class clusterSequence {

  public double sq(double x) { return x*x; }
  
  private boolean kt_alg;
  private double  jetR2;
  private int     num;

  /**
   * pseudoJet class
   */
  private class pseudoJet {
    public double px, py, pz, E, rap, phi, diB;
    public int i;

    public pseudoJet(double px, double py, double pz, double E) {
      this.px = px;
      this.py = py;
      this.pz = pz;
      this.E  = E;
      this.i  = num++;

      rap = 0.5*Math.log((E+pz)/(E-pz));
      phi = (px == 0. && py == 0. ? 0. : Math.atan2(py,px));

      diB = (kt_alg ? pt2() : 1./pt2());
    }

    public pseudoJet(pseudoJet a, pseudoJet b) {
      this( a.px + b.px, a.py + b.py,
            a.pz + b.pz, a.E  + b.E );
    }

    public double pt2() { return px*px + py*py; }
    
    public double dij(pseudoJet p) {
      double deltaPhi = Math.abs(phi-p.phi);
      if (deltaPhi > Math.PI) deltaPhi = 2*Math.PI - deltaPhi;
      return Math.min(diB,p.diB)*( sq(rap-p.rap) + sq(deltaPhi) );
    }

    public String toString() {
      // return "(" + px + ", " + py + ", " + pz + ", " + E + ")";
      return Double.toString(Math.sqrt(pt2())); // print pT
    }
  }

  /**
   * Constructor
   */
  public clusterSequence(boolean kt_alg, double jetR) {
    this.kt_alg = kt_alg;
    this.jetR2  = jetR*jetR;
  }
  
  /**
   * clustering function
   */
  public List<ParticleD> cluster(List<ParticleD> particles) {
    num = 0;
    final int n = particles.size();
    int n_ok = n;
  
    pseudoJet[] pp = new pseudoJet[n];
    boolean[]   ok = new boolean[n];
    double[][] dij = new double[n][n];
    List<ParticleD> jets = new ArrayList<ParticleD>();
  
    // read input particles
    for (int i=0; i<n; ++i) {
      pp[i] = new pseudoJet(
        particles.get(i).px(), particles.get(i).py(),
        particles.get(i).pz(), particles.get(i).e ()
      );
      ok[i] = true;
    }

    // cache pairwise distances
    for (int i=1, k=0; i<n; ++i)
      for (int j=0; j<i; ++j, ++k)
        dij[i][j] = pp[i].dij( pp[j] );
        
    int i1 = 1, i2 = 0;
        
    while (n_ok > 0) {
    
      double  dist  = Double.MAX_VALUE;
      boolean merge = false;

      // find smallest single distance
      for (int i=0; i<n; ++i) {
        if (!ok[i]) continue;
        double d = pp[i].diB*jetR2;
        if (d < dist) {
          dist = d;
          i1 = i;
        }
      }
      
      // find closest pair
      for (int i=1, k=0; i<n; ++i) {
        if (!ok[i]) {
          k += i;
          continue;
        }
        for (int j=0; j<i; ++j, ++k) {
          if (!ok[j]) continue;

          double d = dij[i][j];

          if (d < dist) {
            dist = d;
            i1 = i;
            i2 = j;
            if (!merge) merge = true;
          }
        }
      }

      if (merge) {
        // merge particles
        pseudoJet p = new pseudoJet( pp[i1], pp[i2] );

        // print clustering step
        //System.out.format("%3d: merged %3d & %3d | d = %.8e\n",
        //  p.i, pp[i1].i, pp[i2].i, dist);
        
        // "remove" merge particles
        pp[i1] = p;
        ok[i2] = false;

        // cache new pairwise distances
        for (int i=0; i<i1; ++i) {
          if (!ok[i]) continue;
          dij[i1][i] = pp[i].dij( pp[i1] );
        }
        for (int i=i1+1; i<n; ++i) {
          if (!ok[i]) continue;
          dij[i][i1] = pp[i].dij( pp[i1] );
        }

      } else {
        pseudoJet p = pp[i1];

        // identify as jet
        jets.add( new ParticleD( p.px, p.py, p.pz, p.E ) );
        
        //System.out.format("%3d is a Jet | d = %.8f\n", p.i, dist);

        // "remove"
        ok[i1] = false;
      }

      --n_ok;
    }
    
    return jets;
  }
}
