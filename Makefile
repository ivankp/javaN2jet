CLASS = cluster.class

.PHONY: all clean

all: $(CLASS)

$(CLASS): %.class: %.java
	javac $^

clean:
	rm -f $(wildcard *.class)
