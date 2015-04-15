CLASS = cluster.class

.PHONY: all clean

all: $(CLASS)

cluster.class: ParticleD.java cluster.java clusterSequence.java
	javac $^

clean:
	@rm -fv *.class
