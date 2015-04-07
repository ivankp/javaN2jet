import java.util.*;
import java.io.*;

class particle {
  public double px, py, pz, E, pt2, rap, phi;
  public int id;
  private static int num = 0;

  public particle(double px, double py, double pz, double E) {
    this.px = px;
    this.py = py;
    this.pz = pz;
    this.E  = E;

    pt2 = px*px + py*py;
    rap = 0.5*Math.log((E+pz)/(E-pz));
    phi = (px == 0. && py == 0. ? 0. : Math.atan2(py,px));
    id  = ++num;
  }

  public double diB(boolean kt_alg) { return (kt_alg ? pt2 : 1./pt2); }
  public double dij(particle p, double R, boolean kt_alg) {
    double deltaPhi = Math.abs(phi-p.phi);
    if (deltaPhi > Math.PI) deltaPhi = 2*Math.PI - deltaPhi;
    return Math.min(diB(kt_alg),p.diB(kt_alg))
           * ( Math.pow(rap-p.rap,2) + Math.pow(deltaPhi,2) ) / ( R*R );
  }

  public String toString() {
    return "(" + px + ", " + py + ", " + pz + ", " + E + ")";
  }
}

class cluster {
  private static void usage() {
    System.out.println("Usage: java cluster (anti)kt file");
    System.exit(1);
  }

  public static void main(String[] args) throws IOException {
    if (args.length!=2) usage();

    boolean kt_alg = false;
    if (args[0].equals("kt")) kt_alg = true;
    else if (!args[0].equals("antikt")) usage();

    Vector<particle> jets = new Vector<particle>();

    Scanner sc = new Scanner(new File(args[1]));
    while (sc.hasNextDouble()) {
      jets.add( new particle(
        sc.nextDouble(), sc.nextDouble(), sc.nextDouble(), sc.nextDouble()
      ) );
      System.out.println(jets.lastElement());
    }
  }
}
