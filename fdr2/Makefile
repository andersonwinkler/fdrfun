# Regular use, set compiler name, compiler flags, openmp flags zlib flags
CNAME=gcc
#linker flags
LFLAGS=-lm -lz
#c flags
CFLAGS=-O3 -ffast-math -DHAVE_ZLIB

ifeq ($(CXX),g++)
	CFLAGS += -std=gnu99
endif

all:
	$(CNAME) $(CFLAGS) -o fdr2 main.c core.c bky.c niftilib/nifti2_io.c znzlib/znzlib.c -I./niftilib -I./znzlib $(LFLAGS)
