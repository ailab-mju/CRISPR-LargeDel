# g++-7 was the original compiler; any g++/clang++ produces byte identical
# output (see FIXES.md), so fall back to whatever is installed.
# Override with:  make CPP=clang++   (CPP is a built-in make variable, so it
# has to be assigned with '=' rather than '?=' for this default to take.)
CPP = $(shell command -v g++-7 2>/dev/null || command -v g++ 2>/dev/null)
# struct data collides with std::data under -std=gnu++17 (the default from
# g++ 6 onwards), which is why the language level is pinned here.
STD = -std=gnu++14

# ---------------------------------------------------------------- htslib
# With htslib, kmer_align writes coordinate-sorted BAM + .bai by itself and
# needs no samtools at run time.  Without it the same code still sorts, but
# emits SAM.  Point HTS_CFLAGS / HTS_LDFLAGS at a non-standard prefix:
#   make HTS_CFLAGS=-I$$HOME/tools/include \
#        HTS_LDFLAGS='-L$$HOME/tools/lib -Wl,-rpath,$$HOME/tools/lib'
# (an unpacked .deb has no libhtscodecs.so symlink, so add -lhtscodecs there:
#        make HTS_LIBS='-lhts -lhtscodecs')
# or turn it off entirely with:  make HTS=0
HTS ?= auto
HASH := \#
HTS_CFLAGS ?= $(shell pkg-config --cflags htslib 2>/dev/null)
HTS_LDFLAGS ?=
HTS_LIBS ?= -lhts
ifeq ($(HTS),0)
  HTS_FOUND = 0
else
  HTS_FOUND := $(shell printf '$(HASH)include <htslib/hts.h>\nint main(){return hts_version()==0;}\n' > .hts_probe.cc; \
     $(CPP) $(STD) $(HTS_CFLAGS) .hts_probe.cc -o /dev/null $(HTS_LDFLAGS) $(HTS_LIBS) >/dev/null 2>&1 \
     && echo 1 || echo 0; rm -f .hts_probe.cc)
endif
ifeq ($(HTS_FOUND),1)
  HTS_DEF = -DUSE_HTSLIB
  HTS_LINK = $(HTS_LDFLAGS) $(HTS_LIBS)
else
  HTS_DEF =
  HTS_LINK =
endif

FLAGS = -O3 -g -pthread $(STD) $(HTS_CFLAGS) $(HTS_DEF)
OBJ = align.o align_wrap.o arg_parse.o func.o post.o samout.o report.o kmer_align.o

kmer_align: $(OBJ)
	@if [ "$(HTS_FOUND)" = "1" ]; then echo "-- htslib found: BAM output enabled"; \
	 else echo "-- htslib NOT found: SAM output only (install libhts-dev for BAM)"; fi
	$(CPP) $(FLAGS) -o kmer_align $(OBJ) $(HTS_LINK)

# Every object depends on every header.  The per-object lists below were
# incomplete (align.o and post.o did not list arg_parse.h, although both use
# struct Args and post_process() takes it BY VALUE), so changing Args left
# stale objects behind and the binary corrupted its heap at run time
# ("free(): invalid size").  The project is small enough that rebuilding
# everything on any header change is the right trade.
HDR = $(wildcard *.h) dkm.hpp
$(OBJ): $(HDR) Makefile

kmer_align.o: align_wrap.h report.h kmer_align.cpp
	$(CPP) $(FLAGS) -c kmer_align.cpp
align_wrap.o: func.h align.h arg_parse.h samout.h report.h align_wrap.cpp
	$(CPP) $(FLAGS) -c align_wrap.cpp
align.o: align.h func.h align.cpp
	$(CPP) $(FLAGS) -c align.cpp
arg_parse.o: arg_parse.h samout.h arg_parse.cpp
	$(CPP) $(FLAGS) -c arg_parse.cpp
post.o: post.h samout.h post.cpp
	$(CPP) $(FLAGS) -c post.cpp
func.o: func.h samout.h func.cpp
	$(CPP) $(FLAGS) -c func.cpp
samout.o: samout.h func.h arg_parse.h report.h samout.cpp
	$(CPP) $(FLAGS) -c samout.cpp
report.o: report.h report_tpl.h func.h arg_parse.h report.cpp
	$(CPP) $(FLAGS) -c report.cpp

.PHONY: clean
clean:
	rm -f *.o kmer_align .hts_probe.cc
