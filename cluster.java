import java.util.*;
import java.io.*;

class pseudoJet {
  public static boolean kt_alg;
  public static double jetR;

  public double px, py, pz, E, rap, phi, diB, dij;
  public int i, // own index
             j; // nearest neighbor index

  pseudoJet m1, m2; // mothers

  public pseudoJet(double px, double py, double pz, double E) {
    this.px = px;
    this.py = py;
    this.pz = pz;
    this.E  = E;

    rap = 0.5*Math.log((E+pz)/(E-pz));
    phi = (px == 0. && py == 0. ? 0. : Math.atan2(py,px));

    diB = (kt_alg ? pt2() : 1./pt2());
    dij = Double.MAX_VALUE;
  }

  public pseudoJet(pseudoJet a, pseudoJet b) {
    this( a.px + b.px, a.py + b.py,
          a.pz + b.pz, a.E  + b.E  );

    m1 = a;
    m2 = b;
  }

  public double pt2() { return px*px + py*py; }

  public void update_dij(pseudoJet p) {
    double deltaPhi = Math.abs(phi-p.phi);
    if (deltaPhi > Math.PI) deltaPhi = 2*Math.PI - deltaPhi;
    double dik = Math.min(diB,p.diB)
      * ( Math.pow(rap-p.rap,2) + Math.pow(deltaPhi,2) ) / ( jetR*jetR );

    if (dik < dij) { dij = dik; j = p.i; }
  }

  public String toString() {
    // return "(" + px + ", " + py + ", " + pz + ", " + E + ")";
    return i + ": " + diB;
  }
}

class cluster {
  private static void usage() {
    System.out.println("Usage: java cluster (anti)kt R file");
    System.exit(1);
  }

  public static void main(String[] args) throws IOException {
    if (args.length!=3) usage();

    pseudoJet.kt_alg = false;
    if (args[0].equals("kt")) pseudoJet.kt_alg = true;
    else if (!args[0].equals("antikt")) usage();

    pseudoJet.jetR = Double.parseDouble(args[1]);

    Vector<pseudoJet> jets = new Vector<pseudoJet>();

    Scanner sc = new Scanner(new File(args[2]));
    while (sc.hasNextDouble()) {
      jets.add( new pseudoJet(
        sc.nextDouble(), sc.nextDouble(), sc.nextDouble(), sc.nextDouble()
      ) );
      jets.lastElement().i = jets.size()-1;
      System.out.println(jets.lastElement());
    }
  }
}
