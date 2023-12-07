ifeq ($(shell uname),Linux)
  CC=/usr/bin/gcc
	SHEXT=so
else   # Darwin
  CC=/usr/bin/clang
	SHEXT=dylib
endif


CFLAGS	 := -O3 -Wall -Wextra -Wpedantic -fno-omit-frame-pointer
SRC_DIR = ./src
INC_DIR = ./include

SHARED_FLAGS := -fPIC -shared

INCPATH  := -I$(INC_DIR)
SRCS :=  $(wildcard $(SRC_DIR)/*.c)

main_test: tests/main_test.c $(SRCS) 
	$(CC) $(CFLAGS) $(INCPATH) $(SOURCES) $^ -o $@

bench: tests/bench.c $(SRCS) 
	$(CC) $(CFLAGS) $(INCPATH) $(SOURCES) $^ -o $@	

lib: $(SRCS) 
	$(CC) $(CFLAGS) $(SHARED_FLAGS) $(INCPATH) $(SOURCES) $^ -o libpolyntt.$(SHEXT)

.PHONY: clean

clean:
	-rm main_test bench libpolyntt.$(SHEXT)