CC=gcc
CXX=g++

CXXFLAGS= -std=c++11 -O2
CFLAGS = -O2

.PHONY: all
all: bin/gvcfgenotyper

HTSDIR=external/htslib-1.6
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS = -I$(HTSDIR) -Isrc/cpp/ -Isrc/c/
LFLAGS = -lz -lm -lpthread

debug: CXXFLAGS = -g -O1 -Wall
debug: CFLAGS = -g -O1 
debug: all

profile: CXXFLAGS = -pg -O2 
profile: CFLAGS =  -pg -O2 
profile: all

##generates a version
ifneq "$(wildcard .git)" ""
VERSION = $(shell git describe --always)
endif
build/version.hh:
	echo '#define GG_VERSION "$(VERSION)"' > $@

OBJS=$(shell for i in src/cpp/*.cpp;do echo build/$$(basename $${i%cpp})o;done)
OBJS+=$(shell for i in src/c/*.c;do echo build/$$(basename $${i%c})o;done)

-include $(addsuffix .d,$(OBJS) )

build/%.o: src/cpp/%.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o $@ $<
	$(CXX) -MT $@ -MM $(CXXFLAGS) $(IFLAGS) $< -o $@.d
build/%.o: src/c/%.c 
	$(CC)  $(CFLAGS) $(IFLAGS) -c -o $@ $<
	$(CC) -MT $@ -MM $(CFLAGS) $(IFLAGS) $< -o $@.d

bin/gvcfgenotyper: build/version.hh $(OBJS) $(HTSLIB)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(IFLAGS) $(HTSLIB) $(LFLAGS) $(CXXFLAGS)
.PHONY: clean
clean:
	rm build/* bin/gvcfgenotyper
