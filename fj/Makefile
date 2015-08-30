FJ_DIR    := $(shell fastjet-config --prefix)
FJ_CFLAGS := -I$(FJ_DIR)/include
FJ_LIBS   := -L$(FJ_DIR)/lib -lfastjet

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)

.PHONY: all clean

EXE := cmp one event

all: $(EXE)

one event: %: %.cc
	@echo CC $@
	@g++ -std=c++11 -Wall -O3 $(FJ_CFLAGS) $^ -o $@ $(FJ_LIBS)

cmp: %: %.cc
	@echo CC $@
	@g++ -std=c++11 -Wall -O3 $(FJ_CFLAGS) $(ROOT_CFLAGS) $^ -o $@ $(FJ_LIBS) $(ROOT_LIBS)

clean:
	@rm -fv $(EXE)
