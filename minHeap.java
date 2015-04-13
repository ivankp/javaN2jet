import java.util.*;
import java.io.*;

abstract class heapNode {
  int hi; // heap index
  Double obj;
  abstract double val();
  public heapNode(Double x) { obj = x; }
}

class minHeap {
  private heapNode[] heap;
  private int last;
  
  public minHeap(int size) {
    heap = new heapNode[size];
    last = -1;
  }

/*
  private int child1(int i) { return (i<<1) + 1; }
  private int child2(int i) { return (i<<1) + 2; }
  private int parent(int i) { return ((i-1)>>1); }
*/
  
  public String toString() {
    final int nl = binlog(last+1);
    String out = new String("\n");

    int h = 0;
    for (int l=0; l<=nl; ++l) {
      for (int e=0; e<Math.pow(2,l); ++e) {
        if (e>0) for (int s=0; s<Math.pow(2,nl-l+1)-1; ++s) out += "  ";
        else     for (int s=0; s<Math.pow(2,nl-l  )-1; ++s) out += "  ";
        out += String.format("%2.0f",heap[h++].val());
        if (h>last) break;
      }
      out += "\n";
      if (h>last) break;
    }

    return out;
  }
  
  private void swap(int i, int j) {
    System.out.printf("%2d (%2.0f) <--> %2d (%2.0f)\n",i,heap[i].val(),j,heap[j].val());
  
    // swap nodes
    heapNode a  = heap[i];
    heap[i] = heap[j];
    heap[j] = a;
    
    // swap heap indices
    heap[i].hi = heap[i].hi + heap[j].hi;
    heap[j].hi = heap[i].hi - heap[j].hi;
    heap[i].hi = heap[i].hi - heap[j].hi;
  }
  
  private int sift_up(int i) { // swap i with its parent
    if (i==0) return 0;
    final int p = ((i-1)>>1); // parent

    if (heap[i].val() < heap[p].val()) {
      swap(i,p);
      return p;
    } else return i;
  }
  
  private int sift_down(int i) { // swap i with its child
    if (i==0) return 0;
    final int c1 = (i<<1) + 1, c2 = c1+1; // children
    final int c  = (heap[c1].val() < heap[c2].val() ? c1 : c2);

    if (heap[c].val() < heap[i].val()) {
      swap(i,c);
      return c;
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
  
  public void insert(heapNode a) {
    heap[++last] = a;
    if (last == 0) a.hi = last;
    else {
      int i=last, j=sift_up(i);
      while (i!=j) {
        i = j;
        j = sift_up(i);
      }
    }
  }
  
  public heapNode remove(int i) {
    if (last<i) return null;
    else {
      heapNode ret = heap[i];
      heap[i]  = heap[--last];
      heap[i].hi = i;
      
      int k=i, j=sift_down(k);
      while (k!=j) {
        k = j;
        j = sift_down(k);
      }
      
      return ret;
    }
  }
  
  public heapNode pop() { return remove(0); }
}
