import java.util.*;
import java.io.*;

class clusterSequence {

  public static double sq(double x) { return x*x; }
  
  private boolean kt_alg;
  private double  jetR2;
  private int     num;

  /**
   * pseudoJet class
   */
  private class pseudoJet {
    public double px, py, pz, E, rap, phi, diB, dij;
    public int id, i, j;

    public pseudoJet(double px, double py, double pz, double E, int i) {
      this.px = px;
      this.py = py;
      this.pz = pz;
      this.E  = E;
      this.id = num++;
      this.i  = i;
      this.j  = i;

      rap = 0.5*Math.log((E+pz)/(E-pz));
      phi = (px == 0. && py == 0. ? 0. : Math.atan2(py,px));

      diB = (kt_alg ? pt2() : 1./pt2());
      dij = Double.MAX_VALUE;
    }

    public pseudoJet(pseudoJet a, pseudoJet b, int i) {
      this( a.px + b.px, a.py + b.py,
            a.pz + b.pz, a.E  + b.E, i);
    }

    public double pt2() { return px*px + py*py; }
    
    public double update_dij(pseudoJet p) {
      double deltaPhi = Math.abs(phi-p.phi);
      if (deltaPhi > Math.PI) deltaPhi = 2*Math.PI - deltaPhi;
      double dik = Math.min(diB,p.diB)*( sq(rap-p.rap) + sq(deltaPhi) )/jetR2;
      
      if (dik < dij) { dij = dik; j = p.i; }
      
      return dik;
    }

    /*public String toString() {
      // return "(" + px + ", " + py + ", " + pz + ", " + E + ")";
      return Double.toString(Math.sqrt(pt2())); // print pT
    }*/
  }
    
  private void update_dij(pseudoJet a, pseudoJet b) {
    double deltaPhi = Math.abs(a.phi-b.phi);
    if (deltaPhi > Math.PI) deltaPhi = 2*Math.PI - deltaPhi;
    double dik = Math.min(a.diB,b.diB)*( sq(a.rap-b.rap) + sq(deltaPhi) )/jetR2;
    
    if (dik < a.dij) { a.dij = dik; a.j = b.i; }
    if (dik < b.dij) { b.dij = dik; b.j = a.i; }
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
    // double[][] dij = new double[n][n];
    List<ParticleD> jets = new ArrayList<ParticleD>();
  
    // read input particles
    for (int i=0; i<n; ++i) {
      pp[i] = new pseudoJet(
        particles.get(i).px(), particles.get(i).py(),
        particles.get(i).pz(), particles.get(i).e (), i
      );
    }

    // calculate minimum pairwise distances
    for (int i=1; i<n; ++i)
      for (int j=0; j<i; ++j)
        update_dij(pp[i],pp[j]);
        
    int i1 = 1, i2 = 0;
    boolean merge = false;
    
    // loop until pseudoJet are used up
    while (n_ok > 0) {
    
      double dist = Double.MAX_VALUE;
      
      // find smallest distance
      for (int i=0; i<n; ++i) {
        if (pp[i]!=null) {
          if (pp[i].diB < dist) { dist = pp[i].diB; merge = false; i1 = i; }
          if (pp[i].dij < dist) { dist = pp[i].dij; merge = true;  i1 = i; i2 = pp[i].j; }
        }
      }
      
      //System.out.format("%3d %3d  dist = %.8e\n", i1, i2, dist);
      
      //if (pp[i1]==null || pp[i2]==null) System.exit(0);

      // Either merge or identify a jet
      if (merge) {
      
        // prefer elements with lower indices
        /*if (i1 > i2) { // swap
          i1 = i1 + i2;
          i2 = i1 - i2;
          i1 = i1 - i2;
        }*/
      
        // merge particles
        pseudoJet p = new pseudoJet( pp[i1], pp[i2], i1 );

        // print clustering step
        //System.out.format("%3d: merged %3d & %3d | d = %.8e\n",
        //  p.id, pp[i1].i, pp[i2].i, dist);
        
        // "remove" merge particles
        pp[i1] = p;
        pp[i2] = null;

        // recompute pairwise distances
        for (int i=0; i<n; ++i) {
          if (pp[i]==null) continue;
          if (pp[i].j==i1 || pp[i].j==i2) {
            pp[i].dij = Double.MAX_VALUE;
            for (int j=0; j<n; ++j) {
              if (i==j) continue;
              if (pp[j]==null) continue;
              pp[i].update_dij(pp[j]);
            }
          }
        }

      } else {
        pseudoJet p = pp[i1];

        // identify as jet
        jets.add( new ParticleD( p.px, p.py, p.pz, p.E ) );
        
        // System.out.format("%3d (%3d) is a Jet | d = %.8f\n", p.id, p.i, dist);

        // "remove"
        pp[i1] = null;

        // recompute pairwise distances
        for (int i=0; i<n; ++i) {
          if (pp[i]==null) continue;
          if (pp[i].j==i1) {
            pp[i].dij = Double.MAX_VALUE;
            for (int j=0; j<n; ++j) {
              if (pp[j]==null || j==i) continue;
              pp[i].update_dij(pp[j]);
            }
          }
        }
        
      }
      
      --n_ok;
    }
    
    return jets;
  }
}
