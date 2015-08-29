.PHONY: all clean

all: cluster.class test2.class

ifeq ($(shell \
		echo `which fastjet-config | rev | cut -d'/' -f1 && \
		 			which root | rev | cut -d'/' -f1` | rev \
	), root fastjet-config)

all: fj_cmp

FJ_DIR    := $(shell fastjet-config --prefix)
FJ_CFLAGS := -I$(FJ_DIR)/include
FJ_LIBS   := -L$(FJ_DIR)/lib -lfastjet

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)

endif

%.class: %.java
	javac *.java

fj_cmp: %: %.cc
	g++ -std=c++11 -Wall -O3 $(FJ_CFLAGS) $(ROOT_CFLAGS) $^ -o $@ $(FJ_LIBS) $(ROOT_LIBS)

clean:
	@rm -fv *.class
