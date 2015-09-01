.PHONY: all clean

all: Cluster.class Test.class Benchmark.class

%.class: %.java
	javac *.java

Cluster.class Test.class Benchmark.class: ClusterSequence.java ParticleD.java

clean:
	@rm -fv *.class
