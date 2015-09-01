.PHONY: all clean

all: Cluster.class Test.class

%.class: %.java
	javac *.java

Cluster.class Test.class: ClusterSequence.java ParticleD.java

clean:
	@rm -fv *.class
