CLASS = cluster.class

.PHONY: all clean

all: $(CLASS)

cluster.class: ParticleD.java minHeap.java clusterSequence.java cluster.java
	javac $^

clean:
	rm -f $(wildcard *.class)
