import java.util.*;
import java.io.*;

/**
 * @author Ivan Pogrebnyak
 */

class cluster {
  public static void main(String[] args) throws IOException {
    // check arguments
    if (args.length!=3) {
      System.out.println("Usage: java cluster [kt,antikt,cambridge] R file");
      System.exit(1);
    }

    // set algorithm type and jet radius, R.
    clusterSequence seq = new clusterSequence(
      args[0], Double.parseDouble(args[1])
    );

    List<ParticleD> pp = new ArrayList<ParticleD>();

    // read input file and collect input particles
    Scanner dat = new Scanner(new File(args[2]));
    while (dat.hasNextDouble()) {
      pp.add( new ParticleD(
        dat.nextDouble(), dat.nextDouble(), dat.nextDouble(), dat.nextDouble()
      ) );
      // System.out.println(pp.get(pp.size()-1));
    }

    // perform jet clustering
    List<ParticleD> jets = null;
    for (int i=0; i<10; ++i) {
      long startTime = System.currentTimeMillis();
      jets = seq.cluster(pp);
      System.out.println("N2 run time: " +
        (System.currentTimeMillis()-startTime) + " ms"
      );
    }

    // print
    Collections.sort(jets);
    for (ParticleD j: jets) System.out.format("%.8e\n",j.perp());
  }
}
