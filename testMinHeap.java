import java.util.*;
import java.io.*;

class doubleHeapNode extends heapNode {
  double x;
  doubleHeapNode(double x) { this.x = x; }
  double val() { return x; }
}

class testMinHeap {
  public static void main(String[] args) throws IOException {
    minHeap heap = new minHeap(100);
/*
    heap.insert(new doubleHeapNode(60));
    heap.insert(new doubleHeapNode(14));
    heap.insert(new doubleHeapNode( 5));
    heap.insert(new doubleHeapNode(23));
    heap.insert(new doubleHeapNode( 2));
    heap.insert(new doubleHeapNode(10));
    heap.insert(new doubleHeapNode(15));
*/
    doubleHeapNode a = new doubleHeapNode(10);

    heap.insert(new doubleHeapNode( 1));
    heap.insert(a);
    heap.insert(new doubleHeapNode( 5));
    heap.insert(new doubleHeapNode(20));
    heap.insert(new doubleHeapNode(30));
    heap.insert(new doubleHeapNode( 6));
    
    System.out.println(heap);
    
    heap.remove(3);
    
    System.out.println(heap);
    
    heap.pop();
    
    System.out.println(heap);
    
    a.x = 1;
    
    heap.update(a.hi);
    
    System.out.println(heap);
    
/*
    heap.pop();
    
    System.out.println(heap);
    
    heap.remove(1);
    
    System.out.println(heap);
    
    heap.insert(new doubleHeapNode(25));
    
    System.out.println(heap);
    
    heap.remove(3);
    
    System.out.println(heap);
*/
/*
    String out = new String();
    heapNode d;
    while ((d=heap.pop())!=null) {
      System.out.println(heap);
      out += String.format(" %2.0f",d.val());
    }
    System.out.println(out);
*/
  }
}
