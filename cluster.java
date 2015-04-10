import java.util.*;
import java.io.*;

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
    clusterSequence seq = new clusterSequence(
      kt_alg, Double.parseDouble(args[1])
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
    for (ParticleD j: jets) System.out.format("%.8f\n",j.perp());
  }
}
