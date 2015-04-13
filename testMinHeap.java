import java.util.*;
import java.io.*;

class doubleHeapNode extends heapNode {
  doubleHeapNode(double x) { super(new Double(x)); }
  double val() { return obj.doubleValue(); }
}

class testMinHeap {
  public static void main(String[] args) throws IOException {
    minHeap test = new minHeap(100);
    test.insert(new doubleHeapNode(50));
    test.insert(new doubleHeapNode(60));
    test.insert(new doubleHeapNode(50));
    test.insert(new doubleHeapNode(14));
    test.insert(new doubleHeapNode( 5));
    test.insert(new doubleHeapNode(23));
    test.insert(new doubleHeapNode( 2));
    
    System.out.println(test);
  }
}
