import java.util.*;
import java.io.*;

class node {
  int hi; // heap index
  int data;

  public node(int x) { data = x; }
}

class minHeap {
  private node[] heap;
  private int last;
  
  public minHeap(int size) {
    heap = new node[size];
    last = -1;
  }
  
  private int child1(int i) { return (i<<1) + 1; }
  private int child2(int i) { return (i<<1) + 2; }
  private int parent(int i) { return ((i-1)>>1); }
  
  // swap i into j
  private int trickle(int i) {
    if (i==0) return 0;
    final int p = parent(i);
    if (heap[i].data < heap[p].data) {
      System.out.printf("%2d (%2d) <--> %2d (%2d)\n",i,heap[i].data,p,heap[p].data);
      // swap nodes
      node a  = heap[i];
      heap[i] = heap[p];
      heap[p] = a;
      
      // swap heap indices
      heap[i].hi = heap[i].hi + heap[p].hi;
      heap[p].hi = heap[i].hi - heap[p].hi;
      heap[i].hi = heap[i].hi - heap[p].hi;
      
      return p;
    } else return i;
  }
  
  public static int binlog( int bits ) { // returns 0 for bits=0
    int log = 0;
    if( ( bits & 0xffff0000 ) != 0 ) { bits >>>= 16; log = 16; }
    if( bits >= 256 ) { bits >>>= 8; log += 8; }
    if( bits >= 16  ) { bits >>>= 4; log += 4; }
    if( bits >= 4   ) { bits >>>= 2; log += 2; }
    return log + ( bits >>> 1 );
  }
  
  public void insert(node a) {
    heap[++last] = a;
    if (last == 0) a.hi = last;
    else {
      int i=last, j=trickle(i);
      while (i!=j) {
        i = j;
        j = trickle(i);
      }
    }
    
    for (int i=last, j=-1; i!=j; j=trickle(i)) { }
  }
  
  public String toString() {
    final int nl = binlog(last+1);
    String out = new String("\n");
    //for (int i=0;i<=last;++i)
    //  out += String.format(" %2d",heap[i].data);

    int h = 0;
    for (int l=0; l<=nl; ++l) {
      //out += String.format("%2d   ",l);
      for (int e=0; e<Math.pow(2,l); ++e) {
        if (e>0) for (int s=0; s<Math.pow(2,nl-l+1)-1; ++s) out += "  ";
        else     for (int s=0; s<Math.pow(2,nl-l  )-1; ++s) out += "  ";
        out += String.format("%2d",heap[h++].data);
        if (h>last) break;
      }
      out += "\n";
      if (h>last) break;
    }

    return out;
  }
  
  public static void main(String[] args) throws IOException {
    minHeap test = new minHeap(100);
    test.insert(new node(50));
    test.insert(new node(50));
    test.insert(new node(50));
    test.insert(new node(14));
    test.insert(new node( 5));
    test.insert(new node(23));
    test.insert(new node( 2));
    
    System.out.println(test);
  }
}
