import java.util.*;
import java.io.*;

/**
 * @author Ivan Pogrebnyak
 */

class running_stat {
  private int n;
  private double oldM, newM, oldS, newS;

  public running_stat() { n=0; }

  public void clear() { n=0; }
  public void push(double x) {
    n++;

    // See Knuth TAOCP vol 2, 3rd edition, page 232
    if (n == 1) {
      oldM = newM = x;
      oldS = 0.;
    } else {
      newM = oldM + (x - oldM)/n;
      newS = oldS + (x - oldM)*(x - newM);

      // set up for next iteration
      oldM = newM;
      oldS = newS;
    }
  }
  public int num() { return n; }

  double  mean() { return (n > 0) ? newM : 0.; }
  double   var() { return ( (n > 1) ? newS/(n - 1) : 0. ); }
  double stdev() { return Math.sqrt( var() ); }
};

class Benchmark {
  public static void main(String[] args) throws IOException {
    // check arguments
    if (args.length!=4) {
      System.out.println("Usage: java cluster [kt,antikt,cambridge] R Np Nev");
      System.exit(1);
    }

    // set algorithm type and jet radius, R.
    ClusterSequence seq = new ClusterSequence(
      args[0], Double.parseDouble(args[1])
    );

    int np  = Integer.parseInt(args[2]);
    int nev = Integer.parseInt(args[3]);

    running_stat stat = new running_stat();

    for (int ev=0; ev<nev; ++ev) {
      List<ParticleD> pp = new ArrayList<ParticleD>();

      for (int p=0; p<np; ++p) {
        double m, pt,
               eta = 10.*(Math.acos(1.-2.*Math.random())/Math.PI-0.5),
               phi = 2.*Math.PI*Math.random();
        while (Double.isInfinite(pt = 10.-150.*Math.log1p(-Math.random()))) { }
        while (Double.isInfinite(m  =     -20.*Math.log1p(-Math.random()))) { }

        double px = pt*Math.cos(phi),
               py = pt*Math.sin(phi),
               pz = pt*Math.sinh(eta),
               e  = Math.sqrt(px*px+py*py+pz*pz+m*m);

        pp.add( new ParticleD(px,py,pz,e) );
      }

      long startTime = System.nanoTime();
      List<ParticleD> jets = seq.cluster(pp,0.);
      stat.push((System.nanoTime()-startTime)/1000.);
    }

    System.out.printf("%f ± %f μs\n",stat.mean(),stat.stdev());
  }
}
