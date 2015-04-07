import java.util.*;
import java.io.*;

class particle {
  public double px, py, pz, E;
  public particle(double px, double py, double pz, double E) {
    this.px = px;
    this.py = py;
    this.pz = pz;
    this.E  = E;
  }
  public String toString() {
    return "(" + px + ", " + py + ", " + pz + ", " + E + ")";
  }
}

class cluster {
  public static void main(String[] args) throws IOException {
    Scanner sc = new Scanner(new File(args[0]));
    while (sc.hasNextDouble()) {
      particle p = new particle(
        sc.nextDouble(), sc.nextDouble(), sc.nextDouble(), sc.nextDouble()
      );
      System.out.println(p);
    }
  }
}
