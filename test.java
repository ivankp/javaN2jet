import java.util.*;
import java.io.*;

class test {
  public static void main(String[] args) throws IOException {
    // check arguments
    if (args.length!=10) {
      System.out.println("Usage: java cluster [kt,antikt,cambridge] R px1 py1 pz1 E1 px2 py2 pz2 E2");
      System.exit(1);
    }

    // set algorithm type and jet radius, R.
    clusterSequence seq = new clusterSequence(
      args[0], Double.parseDouble(args[1])
    );

    List<ParticleD> pp = new ArrayList<ParticleD>();

    pp.add( new ParticleD(
      Double.parseDouble(args[2]),
      Double.parseDouble(args[3]),
      Double.parseDouble(args[4]),
      Double.parseDouble(args[5])
    ) );
    pp.add( new ParticleD(
      Double.parseDouble(args[6]),
      Double.parseDouble(args[7]),
      Double.parseDouble(args[8]),
      Double.parseDouble(args[9])
    ) );

    // perform jet clustering
    List<ParticleD> jets = null;
    jets = seq.cluster(pp);

    // print
    Collections.sort(jets);
    for (ParticleD j: jets) System.out.format("%.8e\n",j.perp());
  }
}