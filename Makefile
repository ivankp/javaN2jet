.PHONY: all clean

all: cluster.class test2.class

%.class: %.java
	javac *.java

cluster.class test2.class: clusterSequence.java ParticleD.java

clean:
	@rm -fv *.class
