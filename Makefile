.PHONY: all clean

all: Cluster.class Test.class

%.class: %.java
	javac *.java

cluster.class test2.class: ClusterSequence.java ParticleD.java

clean:
	@rm -fv *.class
