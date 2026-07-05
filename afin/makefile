# afin - Fast-Plast seed-and-extend assembler
#
# Conda-build aware: honors CXX, CXXFLAGS, CPPFLAGS, LDFLAGS and PREFIX/DESTDIR
# supplied by the build environment, while still building standalone from a
# plain `make` in a git clone.

CXX      ?= g++
CXXFLAGS ?= -O2 -g
CXXFLAGS += -std=c++11 -pthread
LDLIBS   += -lz

PREFIX  ?= /usr/local
DESTDIR ?=
BINDIR  := $(DESTDIR)$(PREFIX)/bin

TARGET := afin
OBJS := afin.o contig.o process.o read.o revcomp.o mismatch.o log.o gzip.o \
        fusion.o match.o readlist.o extension.o contiglist.o

.PHONY: all clean install
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

install: $(TARGET)
	mkdir -p $(BINDIR)
	install -m 0755 $(TARGET) $(BINDIR)/$(TARGET)

clean:
	rm -f *.o $(TARGET) test
