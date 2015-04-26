import java.util.*;
import java.io.*;

class cluster {
  private static void usage() {
    System.out.println("Usage: java cluster [kt,antikt,cambridge] R file");
    System.exit(1);
  }

  public static void main(String[] args) throws IOException {
    // check arguments
    if (args.length!=3) usage();

    // set algorithm type and jet radius, R.
    jetAlg jet_alg = null;
    if (args[0].equalsIgnoreCase("kt")) jet_alg = jetAlg.kt;
    else if (args[0].equalsIgnoreCase("antikt")) jet_alg = jetAlg.antikt;
    else if (args[0].equalsIgnoreCase("cambridge")) jet_alg = jetAlg.cambridge;
    else {
      System.out.println("Unrecognized clustering algorithm: "+args[0]);
      usage();
    }

    // set up clustering algorithm
    clusterSequence seq = new clusterSequence(
      jet_alg, Double.parseDouble(args[1])
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
      System.out.println("Run time: " +
        (System.currentTimeMillis()-startTime) + " ms"
      );
    }

    // print
    Collections.sort(jets);
    for (ParticleD j: jets) System.out.format("%.8e\n",j.perp());
  }
}
