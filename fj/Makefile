FJ_DIR    := $(shell fastjet-config --prefix)
FJ_CFLAGS := -I$(FJ_DIR)/include
FJ_LIBS   := -L$(FJ_DIR)/lib -lfastjet

.PHONY: all clean

GREP_EXE := grep -r '^ *int *main *(' . | cut -d':' -f1
EXE := $(patsubst ./%.cc,%,$(shell $(GREP_EXE)))

all: $(EXE)

$(EXE): %: %.cc
	@echo CC $@
	@g++ -std=c++11 -Wall -O3 $(FJ_CFLAGS) $^ -o $@ $(FJ_LIBS) -pthread

clean:
	@rm -fv $(EXE)

