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
    public double px, py, pz, E, rap, phi, Rij, diB, dij;
    public int id;
    public pseudoJet prev, next, near;

    public pseudoJet(double px, double py, double pz, double E) {
      this.px = px;
      this.py = py;
      this.pz = pz;
      this.E  = E;
      this.id = num++;

      rap = 0.5*Math.log((E+pz)/(E-pz));
      phi = (px == 0. && py == 0. ? 0. : Math.atan2(py,px));

      diB = (kt_alg ? pt2() : 1./pt2());
      Rij = Double.MAX_VALUE;

      // System.out.printf("%3d: diB %.8e\n",id,diB);
    }

    public void remove() {
      // System.out.printf("removing %3d\n",id);

      if (prev==null) first = next;
      else prev.next = next;

      if (next!=null) next.prev = prev;
    }

    public void merge() {
      first.prev = new pseudoJet(
        px + near.px, py + near.py,
        pz + near.pz, E  + near.E );
      first.prev.next = first;
      // first.prev.near = this; // trick
      first = first.prev;

      this.remove();
      near.remove();
    }

    public double pt2() { return px*px + py*py; }

    public double update_near(pseudoJet p) {
      double deltaPhi = Math.abs(phi-p.phi);
      if (deltaPhi > Math.PI) deltaPhi = 2*Math.PI - deltaPhi;
      double Ril = sq(rap-p.rap) + sq(deltaPhi);

      if (Ril < Rij) { Rij = Ril; near = p; }
      return Ril;
    }

    public double update_near_both(pseudoJet p) {
      double Ril = update_near(p);
      if (Ril < p.Rij) { p.Rij = Ril; p.near = this; }
      return Ril;
    }

    public void update_dij() {
      dij = Math.min(diB,near.diB)*Rij/jetR2;
    }
  }

  private pseudoJet first;

  /**
   * Constructor
   */
  public clusterSequence(boolean kt_alg, double jetR) {
    this.kt_alg = kt_alg;
    this.jetR2  = jetR*jetR;
    first = null;
  }

  /**
   * clustering function
   */
  public List<ParticleD> cluster(List<ParticleD> particles) {
    final int n = particles.size();
    num = 0;

    List<ParticleD> jets = new ArrayList<ParticleD>();
    pseudoJet p;

    if (n==0) return jets;

    // read input particles
    first = new pseudoJet(
      particles.get(0).px(), particles.get(0).py(),
      particles.get(0).pz(), particles.get(0).e ()
    );
    p = first;
    for (int i=1; i<n; ++i) {
      p.next = new pseudoJet(
        particles.get(i).px(), particles.get(i).py(),
        particles.get(i).pz(), particles.get(i).e ()
      );
      p.next.prev = p;
      p = p.next;
    }

    // find nearest neighbors
    for (p=first.next; p!=null; p=p.next)
      for (pseudoJet q=first; q!=p; q=q.next)
        p.update_near_both(q);

    // calculate minimum pairwise distances
    for (p=first; p!=null; p=p.next) {
      p.update_dij();
      // System.out.println(p.dij);
    }

    boolean merge = false;

    // loop until pseudoJet are used up
    while (first != null) {

      double dist = Double.MAX_VALUE;

      // find smallest distance
      for (pseudoJet q=first; q!=null; q=q.next) {
        if (q.diB < dist) { p = q; dist = q.diB; merge = false; }
        if (q.dij < dist) { p = q; dist = q.dij; merge = true;  }
        // System.out.format("%3d: %.8e %3d %.8e\n",q.id,q.diB,q.near.id,q.dij);
      }

      // Either merge or identify a jet
      if (merge) {

        // merge particles
        p.merge();

        // print clustering step
        // System.out.format("%3d & %3d | d = %.5e\n",p.id, p.near.id, dist);
        
        // recompute pairwise distances for the new pseudoJet
        for (pseudoJet q=first.next; q!=null; q=q.next) {
          // the new particle is first
          first.update_near_both(q);
          if (q.near==first) q.update_dij();
        }
        first.update_dij();

        // recompute pairwise distances for the rest
        for (pseudoJet p1=first.next; p1!=null; p1=p1.next) {
          if (p1.near==p || p1.near==p.near) {
            p1.Rij = Double.MAX_VALUE;
            for (pseudoJet p2=first; p2!=null; p2=p2.next) {
              if (p1!=p2) p1.update_near(p2);
            }
            p1.update_dij();
          }
        }

      } else {
        // identify as jet
        jets.add( new ParticleD( p.px, p.py, p.pz, p.E ) );

        // print clustering step
        // System.out.format("%3d Jet   | d = %.5e\n", p.id, dist);

        // "remove"
        p.remove();
        // System.out.printf("p.id = %3d\n",p.id);

        // recompute pairwise distances
        for (pseudoJet p1=first; p1!=null; p1=p1.next) {
          // System.out.printf("%3d near %3d\n",p1.id,p1.near.id);
          if (p1.near==p) {
            p1.Rij = Double.MAX_VALUE;
            for (pseudoJet p2=first; p2!=null; p2=p2.next) {
              // System.out.print(".");
              if (p1!=p2) p1.update_near(p2);
            }
            p1.update_dij();
            // System.out.println();
          }
        }

      }

    }

    return jets;
  }
}
