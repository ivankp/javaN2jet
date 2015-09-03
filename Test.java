import java.util.*;
import java.io.*;

/**
 * @author Ivan Pogrebnyak
 */

class Test {
  public static void main(String[] args) throws IOException {
    // check arguments
    if (args.length<10) {
      System.out.println(
        "Usage: java test2 [kt,antikt,cambridge] " +
        "R px1 py1 pz1 E1 px2 py2 pz2 E2 ...");
      System.exit(1);
    }

    // set algorithm type and jet radius, R.
    ClusterSequence seq = new ClusterSequence(
      args[0], Double.parseDouble(args[1])
    );

    List<ParticleD> pp = new ArrayList<ParticleD>();

    for (int i=0, n=(args.length-2)/4; i<n; ++i) {
      pp.add( new ParticleD(
        Double.parseDouble(args[2+i*4]),
        Double.parseDouble(args[3+i*4]),
        Double.parseDouble(args[4+i*4]),
        Double.parseDouble(args[5+i*4])
      ) );
    }

    // perform jet clustering
    List<ParticleD> jets = seq.cluster(pp,0.);

    // print
    Collections.sort(jets);
    for (ParticleD j: jets) System.out.format("%.8e ",j.perp());
  }

  public static double[] cluster(String alg, double jetR, double[] p) {
    // set algorithm type and jet radius, R.
    ClusterSequence seq = new ClusterSequence(alg,jetR);
    List<ParticleD> pp = new ArrayList<ParticleD>();

    for (int i=0, n=p.length/4; i<n; ++i)
      pp.add( new ParticleD(p[i*4], p[i*4+1], p[i*4+2], p[i*4+3]) );

    // perform jet clustering
    List<ParticleD> jets = seq.cluster(pp,0.);
    Collections.sort(jets);
    double[] pts = new double[jets.size()];
    for (int i=0; i<jets.size(); ++i) pts[i] = jets.get(i).perp();
    return pts;
  }
}
