.PHONY: all clean

all: cluster.class

cluster.class: cluster.java clusterSequence.java ParticleD.java
	javac $^

clean:
	@rm -fv *.class
