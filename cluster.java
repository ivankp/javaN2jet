import java.util.*;
import java.io.*;

class pseudoJet {
  public double px, py, pz, E, rap, phi, diB, dij;
  public int i, // own index
             j; // nearest neighbor index
             
  public pseudoJet(double px, double py, double pz, double E, boolean kt_alg, int i) {
    this.px = px;
    this.py = py;
    this.pz = pz;
    this.E  = E;
    this.i  = i;
    this.j  = i; // trick to evoid extra checks in updating dij

    rap = 0.5*Math.log((E+pz)/(E-pz));
    phi = (px == 0. && py == 0. ? 0. : Math.atan2(py,px));

    diB = (kt_alg ? pt2() : 1./pt2());
    dij = Double.MAX_VALUE;
  }

  public pseudoJet(pseudoJet a, pseudoJet b, boolean kt_alg, int i) {
    this( a.px + b.px, a.py + b.py,
          a.pz + b.pz, a.E  + b.E, kt_alg, i);
  }

  public double pt2() { return px*px + py*py; }
  
  private double sq(double x) { return x*x; }
  
  public void update_dij(pseudoJet p, double jetR2) {
    double deltaPhi = Math.abs(phi-p.phi);
    if (deltaPhi > Math.PI) deltaPhi = 2*Math.PI - deltaPhi;
    double dik = Math.min(diB,p.diB)*( sq(rap-p.rap) + sq(deltaPhi) )/jetR2;

    if (dik < dij) { dij = dik; j = p.i; }
  }

  public String toString() {
    // return "(" + px + ", " + py + ", " + pz + ", " + E + ")";
    return Double.toString(Math.sqrt(pt2())); // print pT
  }
}

class clusterSequence {
  private boolean kt_alg;
  private double  jetR2;

  private List<pseudoJet> pp, jets;
  
  public clusterSequence(boolean kt_alg, double jetR) {
    this.kt_alg = kt_alg;
    this.jetR2  = jetR*jetR;
    
    pp   = new LinkedList<pseudoJet>();
    jets = new LinkedList<pseudoJet>();
  }
  
  public List<pseudoJet> getJets() { return jets; }
  
  public void add(double px, double py, double pz, double E) {
    pseudoJet p = new pseudoJet(px, py, pz, E, kt_alg, pp.size());
    pp.add(p);
    System.out.println(p.i+": "+p);
  }

  public void cluster() {
    System.out.println("\nClustering " + pp.size() + " particles");
  
    for (pseudoJet a: pp)
      for (pseudoJet b: pp)
        if (a!=b) a.update_dij(b,jetR2);
        
    for (pseudoJet a: pp) System.out.print(a.i); System.out.println();
    for (pseudoJet a: pp) System.out.print(a.j); System.out.println();

    int pp_size;

    while ((pp_size = pp.size()) > 0) {

      System.out.println(pp_size);
    
      pseudoJet min_p = null;
      double    min_d = Double.MAX_VALUE;
      boolean   merge = false;

      if (pp_size > 1) {
        for (pseudoJet a: pp) {
          if (a.diB < min_d) { min_d = a.diB; merge = false; min_p = a; }
          if (a.dij < min_d) { min_d = a.dij; merge = true;  min_p = a; }
        }
      
        if (merge) {

          int keep, pop;
          if (min_p.i < min_p.j) { keep = min_p.i; pop = min_p.j; }
          else                   { keep = min_p.j; pop = min_p.i; }

          System.out.println("Merging " + keep + " and " + pop);

          // merge pair
          pp.set(keep, new pseudoJet(min_p, pp.remove(pop), kt_alg, keep));
          
          for (pseudoJet a: pp) System.out.print(a.i); System.out.println();
          for (pseudoJet a: pp) System.out.print(a.j); System.out.println();

          // shift indices because jth element was removed
          for (pseudoJet a: pp) {
            if (a.i > min_p.j) --a.i;
          }

          // recompute pairwise distance
          for (pseudoJet a: pp) {
            if (a.j==min_p.i || a.j==min_p.j) {
              for (pseudoJet b: pp) if (a!=b) {
                a.dij = Double.MAX_VALUE;
                a.update_dij(b,jetR2);
              }
            } else {
              if (a.j > min_p.j) --a.j;
            }
          }
          
          for (pseudoJet a: pp) System.out.print(a.i); System.out.println();
          for (pseudoJet a: pp) System.out.print(a.j); System.out.println();
          
        } else {

          System.out.println("Jet " + min_p.i);
        
          jets.add(min_p);
          pp.remove(min_p.i);
          
          for (pseudoJet a: pp) System.out.print(a.i); System.out.println();
          for (pseudoJet a: pp) System.out.print(a.j); System.out.println();

          // shift indices because ith element was removed
          for (pseudoJet a: pp) {
            if (a.i > min_p.i) --a.i;
          }
          
          // recompute pairwise distance
          for (pseudoJet a: pp) {
            if (a.j==min_p.i) {
              for (pseudoJet b: pp) if (a!=b) {
                a.dij = Double.MAX_VALUE; 
                a.update_dij(b,jetR2);
              }
            } else {
              if (a.j > min_p.i) --a.j;
            }
          }
          
          for (pseudoJet a: pp) System.out.print(a.i); System.out.println();
          for (pseudoJet a: pp) System.out.print(a.j); System.out.println();

        }
        
      } else jets.add(pp.remove(0)); // last pseudoJet is a jet
      
    }
  }
}

class cluster {
  private static void usage() {
    System.out.println("Usage: java cluster (anti)kt R file");
    System.exit(1);
  }

  public static void main(String[] args) throws IOException {
    // check arguments
    if (args.length!=3) usage();

    // set algorithm type, (anti)kt, and jet radius, R.
    boolean kt_alg = false;
    if (args[0].equals("kt")) kt_alg = true;
    else if (!args[0].equals("antikt")) usage();

    // set up clustering algorithm
    clusterSequence jets = new clusterSequence(
      kt_alg, Double.parseDouble(args[1])
    );

    // read input file and collect input particles
    Scanner dat = new Scanner(new File(args[2]));
    while (dat.hasNextDouble()) {
      jets.add(
        dat.nextDouble(), dat.nextDouble(), dat.nextDouble(), dat.nextDouble()
      );
    }

    // perform jet clustering
    jets.cluster();

    // print
    for (pseudoJet j: jets.getJets())
      System.out.println(j);
  }
}
