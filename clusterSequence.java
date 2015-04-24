import java.util.*;
import java.io.*;

class clusterSequence {

  private static double sq(double x) { return x*x; }
  
  private static int mod(int i, int n) {
    return ( (i%=n) < 0 ? i+n : i );
  }
  
  private final static double twopi = 2*Math.PI;

  private boolean kt_alg;
  private double  jetR, jetR2;
  private int     num;

  private pseudoJet first;
  private tileGrid grid;
  
  private int testNumJets() {
    int n = 0;
    for (pseudoJet p=first; p!=null; p=p.next) ++n;
    return n;
  }
  
  // ****************************************************************
  // pseudoJet class ************************************************
  private class pseudoJet {
    public double px, py, pz, E, rap, phi, Rij, diB, dij, RiC1;
    public int id;
    public pseudoJet prev, next, near;
    public tile t;
    public pseudoJet tprev, tnext;

    public pseudoJet(double px, double py, double pz, double E) {
      this.px = px;
      this.py = py;
      this.pz = pz;
      this.E  = E;
      this.id = num++;

      if (E==pz) rap = Double.MAX_VALUE;
      else if (E==-pz) rap = -Double.MAX_VALUE;
      else rap = 0.5*Math.log((E+pz)/(E-pz));
      phi = (px == 0. && py == 0. ? 0. : Math.atan2(py,px)) + Math.PI;

      diB = (kt_alg ? pt2() : 1./pt2());
      Rij = Double.MAX_VALUE;

      // System.out.printf("%3d: diB %.8e\n",id,diB);
    }

    public void remove() {
      // System.out.printf("removing %3d\n",id);

      if (prev==null) first = next;
      else prev.next = next;
      if (next!=null) next.prev = prev;

      if (grid!=null) {
        if (tprev==null) t.first = tnext;
        else tprev.tnext = tnext;
        if (tnext!=null) tnext.tprev = tprev;
      }
    }

    public void merge() {
      first.prev = new pseudoJet(
        px + near.px, py + near.py,
        pz + near.pz, E  + near.E );
      first.prev.next = first;
      first = first.prev;
      
      this.remove();
      near.remove();
    }

    public double pt2() { return px*px + py*py; }

    public boolean update_near(pseudoJet p, boolean both) {
      double deltaPhi = Math.abs(phi-p.phi);
      if (deltaPhi > Math.PI) deltaPhi = twopi - deltaPhi;
      double Ril = sq(rap-p.rap) + sq(deltaPhi);

      if (Ril < Rij) { Rij = Ril; near = p; }
      if (both) {
        if (Ril < p.Rij) { p.Rij = Ril; p.near = this; return true; }
        else return false;
      } else return false;
    }

    public void update_dij() {
      dij = Math.min(diB,near.diB)*Rij/jetR2;
    }
  }
  
  // ****************************************************************
  // tile class *****************************************************
  class tile {
    int irap, iphi;
    double crap, cphi;
    pseudoJet first;
    tile(int irap, int iphi, double crap, double cphi) {
      this.irap = irap;
      this.iphi = iphi;
      this.crap = crap;
      this.cphi = cphi;
    }
  }

  // ****************************************************************
  // tileGrid class *************************************************
  private class tileGrid {
    private final tile[][] tiles;
    private final int nrap, nphi;
    private final double d, r;
    private final double max_rap;
    private final int ng;

    public tileGrid(double R, double max_rap) {
      int _nphi = (int)(twopi/R);
      nphi = (_nphi%2==0 ? _nphi+1 : _nphi );
      ng = (nphi>>1);
      d = twopi/nphi;
      r = d/2;
      int half_nrap = (int)(max_rap/r) + 1;
      this.max_rap = half_nrap*d;
      nrap = 2*half_nrap;
      
      tiles = new tile[nphi][nrap];
      for (int irap=0; irap<nrap; ++irap)
        for (int iphi=0; iphi<nphi; ++iphi)
          tiles[iphi][irap] = new tile(
            irap, iphi,
            (irap-nrap/2)*d + r, iphi*d + r
          );
    }

    public void add(pseudoJet p) {
      int irap = (int)((p.rap+max_rap)/d);
      if (irap < 0) irap = 0;
      else if (irap >= nrap) irap = nrap-1;

      p.t = tiles[ (int)(p.phi/d) ][ irap ];
      
      if (p.t.first==null) {
        p.t.first = p;
      } else {
        p.tnext = p.t.first;
        p.t.first.tprev = p;
        p.t.first = p;
      }

      p.RiC1 = Math.sqrt( sq(p.t.crap-p.rap) + sq(p.t.cphi-p.phi) );
    }

    private void within_tile(pseudoJet p, tile t, boolean both) {
      //System.out.printf("%2d p(%2d,%2d) t(%2d,%2d)\n",p.id,p.t.iphi,p.t.irap,t.iphi,t.irap);
    
      if (t.first!=null) {
        double dphi = 0, drap = 0;
        
        if (t.iphi != p.t.iphi) {
          dphi = Math.abs(p.phi - t.cphi);
          if (dphi > Math.PI) dphi = twopi - dphi;
          dphi -= r;
        }
        if (t.irap != p.t.irap) {
          drap = Math.abs(p.rap - t.crap) - r;
        }
        
        if ( p.Rij >= sq(dphi) + sq(drap) )
          for (pseudoJet q=t.first; q!=null; q=q.tnext)
            if (p.update_near(q, both)) q.update_dij();
      }
    }

    private void within_tile_same(pseudoJet p, boolean both) {
      for (pseudoJet q=p.t.first; q!=null; q=q.tnext)
        if (q!=p)
          if (p.update_near(q, both)) q.update_dij();
    }

    public void update_near(pseudoJet p, boolean both) {
      within_tile_same(p, both);
      if ( p.Rij < sq(p.RiC1) ) return;
    
      for (int k=1; k<=ng; ++k) {
        final int iphi_min = p.t.iphi - k,
                  iphi_max = p.t.iphi + k,
                  irap_min = p.t.irap - k,
                  irap_max = p.t.irap + k;

        int i = iphi_min;

        for (int j=Math.max(irap_min,0); j<=Math.min(irap_max,nrap-1); ++j)
          within_tile(p, tiles[mod(i,nphi)][j], both);

        for (++i; i<iphi_max; ++i) {
          if (irap_min>=0)   within_tile(p, tiles[mod(i,nphi)][irap_min], both);
          if (irap_max<nrap) within_tile(p, tiles[mod(i,nphi)][irap_max], both);
        }

        //if (k>0)
          for (int j=Math.max(irap_min,0); j<=Math.min(irap_max,nrap-1); ++j)
            within_tile(p, tiles[mod(i,nphi)][j], both);

        if ( Math.sqrt(p.Rij) < (p.RiC1 + d*k) ) {
          //System.out.printf("%2d %5d\n",k,testNumJets());
          return;
        }
      }
      
      // left loop
      for (int j=0; j<p.t.irap-ng; ++j)
        for (int i=0; i<nphi; ++i)
          within_tile(p, tiles[i][j], both);
          
      // right loop
      for (int j=p.t.irap+ng+1; j<nrap; ++j)
        for (int i=0; i<nphi; ++i)
          within_tile(p, tiles[i][j], both);

      //System.out.printf("%2d %5d\n",nk,testNumJets());
    }
  }
  
  // ****************************************************************
  // Constructor ****************************************************
  public clusterSequence(boolean kt_alg, double jetR) {
    this.kt_alg = kt_alg;
    this.jetR   = jetR;
    this.jetR2  = jetR*jetR;
    first = null;
  }

  // ****************************************************************
  // clustering function ********************************************
  public List<ParticleD> cluster(List<ParticleD> particles) {
    int n = particles.size();
    num = 0; // start assigning pseudoJet id from 0
    
    // initialize the grid
    if (n>50) grid = new tileGrid(jetR,5);

    List<ParticleD> jets = new ArrayList<ParticleD>();
    pseudoJet p;

    if (n==0) return jets;

    // read input particles -------------------------------
    first = new pseudoJet(
      particles.get(0).px(), particles.get(0).py(),
      particles.get(0).pz(), particles.get(0).e ()
    );
    p = first;
    if (grid!=null) grid.add(p);
    for (int i=1; i<n; ++i) {
      p.next = new pseudoJet(
        particles.get(i).px(), particles.get(i).py(),
        particles.get(i).pz(), particles.get(i).e ()
      );
      p.next.prev = p;
      p = p.next;
      if (grid!=null) grid.add(p);
    }

    // find original nearest neighbors --------------------
    if (grid==null) { // no grid

      for (p=first.next; p!=null; p=p.next)
        for (pseudoJet q=first; q!=p; q=q.next)
          p.update_near(q,true);

    } else { // using grid

      for (p=first; p!=null; p=p.next)
        grid.update_near(p,false);

    }

    // calculate minimum pairwise distances ---------------
    for (p=first; p!=null; p=p.next) {
      p.update_dij();
      // System.out.println(p.dij);
    }

    boolean merge = false;

    // loop until pseudoJets are used up ------------------
    while (first != null) {
    
      if (n<50) grid = null;

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
        
        // the new particle is first
        if (grid!=null) grid.add(first);

        // print clustering step
        // System.out.format("%3d & %3d | d = %.5e\n",p.id, p.near.id, dist);
        
        // recompute pairwise distances
        if (grid==null) { // no grid
        
          // for the new pseudoJet
          for (pseudoJet q=first.next; q!=null; q=q.next)
            if ( first.update_near(q,true) ) q.update_dij();
          first.update_dij();

          // for the rest
          for (pseudoJet p1=first.next; p1!=null; p1=p1.next) {
            if (p1.near==p || p1.near==p.near) {
              p1.Rij = Double.MAX_VALUE;
              for (pseudoJet p2=first; p2!=null; p2=p2.next) {
                if (p1!=p2) p1.update_near(p2,false);
              }
              p1.update_dij();
            }
          }

        } else { // using grid

          // for the new pseudoJet
          grid.update_near(first,true);
          first.update_dij();
          
          // for the rest
          for (pseudoJet p1=first.next; p1!=null; p1=p1.next) {
            if (p1.near==p || p1.near==p.near) {
              p1.Rij = Double.MAX_VALUE;
              grid.update_near(p1,false);
              p1.update_dij();
            }
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
        if (grid==null) { // no grid
        
          for (pseudoJet p1=first; p1!=null; p1=p1.next) {
            // System.out.printf("%3d near %3d\n",p1.id,p1.near.id);
            if (p1.near==p) {
              p1.Rij = Double.MAX_VALUE;
              for (pseudoJet p2=first; p2!=null; p2=p2.next) {
                // System.out.print(".");
                if (p1!=p2) p1.update_near(p2,false);
              }
              p1.update_dij();
              // System.out.println();
            }
          }
          
        } else { // using grid
        
          for (pseudoJet p1=first; p1!=null; p1=p1.next) {
            if (p1.near==p) {
              p1.Rij = Double.MAX_VALUE;
              grid.update_near(p1,false);
              p1.update_dij();
            }
          }
        
        }

      }
      
      --n;

    }

    return jets;
  }
}
