.PHONY: all clean

all: cluster.class

cluster.class: cluster.java clusterSequence.java jetAlg.java ParticleD.java
	javac $^

clean:
	@rm -fv *.class
