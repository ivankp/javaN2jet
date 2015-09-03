FJ_DIR    := $(shell fastjet-config --prefix)
FJ_CFLAGS := -I$(FJ_DIR)/include
FJ_LIBS   := -L$(FJ_DIR)/lib -lfastjet

JDK_INCDIR := $(shell readlink -f `which javac` | sed 's|/bin/javac|/include|')
JRE_LIBDIR := $(shell readlink -f `which java` | sed 's|bin/java|lib/amd64/server|')

.PHONY: all clean

GREP_EXE := grep -r '^ *int *main *(' . | cut -d':' -f1
EXE := $(patsubst ./%.cc,%,$(shell $(GREP_EXE)))

all: $(EXE)

%: %.cc
	@echo CC $@
	@g++ -std=c++11 -Wall -O3 $(FJ_CFLAGS) $^ -o $@ $(FJ_LIBS)

cmp: %: %.cc
	@echo CC $@
	@g++ -std=c++11 -Wall -O3 $(FJ_CFLAGS) $^ -o $@ $(FJ_LIBS) -pthread

cmpjni: %: %.cc
	@echo CC $@
	@g++ -std=c++11 -Wall -g -I$(JDK_INCDIR) $(FJ_CFLAGS) $^ -o $@ \
		-Wl,-rpath=$(JRE_LIBDIR) $(FJ_LIBS) -L$(JRE_LIBDIR) -ljvm

clean:
	@rm -fv $(EXE)
