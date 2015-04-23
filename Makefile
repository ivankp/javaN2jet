CLASS = cluster.class

.PHONY: all clean

all: $(CLASS)

cluster.class: cluster.java clusterSequence.java ParticleD.java
	javac $^

clean:
	@rm -fv *.class
