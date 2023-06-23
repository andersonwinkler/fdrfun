# Regular use, set compiler name, compiler flags, openmp flags zlib flags
CNAME=gcc
#linker flags
LFLAGS=-lm
#c flags
CFLAGS=-O3

ifeq ($(CXX),g++)
	CFLAGS += -std=gnu99
endif

all:
	$(CNAME) $(CFLAGS) -o bky main.c bky.c  $(LFLAGS)
