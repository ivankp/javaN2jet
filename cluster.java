import java.util.*;
import java.io.*;

class pseudoJet {
  public double px, py, pz, E, rap, phi, diB, dij;

  pseudoJet pj; // nearest neighbor
             
  public pseudoJet(double px, double py, double pz, double E, boolean kt_alg) {
    this.px = px;
    this.py = py;
    this.pz = pz;
    this.E  = E;

    rap = 0.5*Math.log((E+pz)/(E-pz));
    phi = (px == 0. && py == 0. ? 0. : Math.atan2(py,px));

    diB = (kt_alg ? pt2() : 1./pt2());
    dij = Double.MAX_VALUE;
  }

  public pseudoJet(pseudoJet a, pseudoJet b, boolean kt_alg) {
    this( a.px + b.px, a.py + b.py,
          a.pz + b.pz, a.E  + b.E, kt_alg );
  }

  public double pt2() { return px*px + py*py; }
  
  private double sq(double x) { return x*x; }
  
  public void update_dij(pseudoJet p, double jetR2) {
    double deltaPhi = Math.abs(phi-p.phi);
    if (deltaPhi > Math.PI) deltaPhi = 2*Math.PI - deltaPhi;
    double dik = Math.min(diB,p.diB)*( sq(rap-p.rap) + sq(deltaPhi) )/jetR2;

    if (dik < dij) { dij = dik; pj = p; }
  }

  public String toString() {
    // return "(" + px + ", " + py + ", " + pz + ", " + E + ")";
    return Double.toString(diB); // print pT
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
    pseudoJet p = new pseudoJet(px, py, pz, E, kt_alg);
    pp.add(p);
    System.out.println(p);
  }

  public void cluster() {
    System.out.println("Clustering " + pp.size() + " particles");
  
    for (pseudoJet a: pp)
      for (pseudoJet b: pp)
        if (a!=b) a.update_dij(b,jetR2);

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
      }
      
      if (merge) {
      
        // merge pair
        pp.add(new pseudoJet(min_p, min_p.pj, kt_alg));
        
        // remove merged particles
        if (!pp.remove(min_p)) System.out.println("not found");
        if (!pp.remove(min_p.pj)) System.out.println("not found");
        
        // recompute pairwise distance
        for (pseudoJet a: pp)
          if (a.pj==min_p || a.pj==min_p.pj || a.pj==null)
            for (pseudoJet b: pp)
              if (a!=b) a.update_dij(b,jetR2);
        
      } else {
      
        jets.add(min_p);
        pp.remove(min_p);
        
        // recompute pairwise distance
        for (pseudoJet a: pp)
          if (a.pj==min_p)
            for (pseudoJet b: pp)
              if (a!=b) a.update_dij(b,jetR2);

      }
    }
  }
}

class cluster {
  Vector<pseudoJet> jets = new Vector<pseudoJet>();

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
