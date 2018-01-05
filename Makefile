.PHONY: all
all: bin/gvcfgenotyper bin/test_gvcfgenotyper

CC=gcc
CXX=g++

CXXFLAGS= -std=c++11 -O2
CFLAGS = -O2

IFLAGS = -Isrc/cpp/lib/ -Isrc/c/
LFLAGS = -lz -lm -lpthread

#testing stuff
TESTFLAGS = -I./external/googletest-release-1.8.0//googletest/include/
include external/googletest-release-1.8.0//googletest/make/Makefile

#htslib stuff
HTSDIR=external/htslib-1.6
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS += -I$(HTSDIR)

#special builds for debugging and profiling
debug: CXXFLAGS = -std=c++11 -g -O1 -Wall
debug: CFLAGS = -g -O1 
debug: all

profile: CXXFLAGS = -std=c++11 -pg
profile: CFLAGS =  -pg
profile: all

##generates a version
ifneq "$(wildcard .git)" ""
VERSION = $(shell git describe --always)
endif
build/version.hh: 
	echo '#define GG_VERSION "$(VERSION)"' > $@

OBJS=$(shell for i in src/cpp/lib/*.cpp;do echo build/$$(basename $${i%cpp})o;done)
OBJS+=$(shell for i in src/c/*.c;do echo build/$$(basename $${i%c})o;done)
TESTOBJS=$(shell for i in src/cpp/test/*.cpp;do echo build/$$(basename $${i%cpp})o;done)

-include $(addsuffix .d,$(OBJS) )
-include $(addsuffix .d,$(TESTOBJS) )

build/%.o: src/cpp/test/%.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) $(TESTFLAGS) -c -o $@ $<
	$(CXX) -MT $@ -MM $(CXXFLAGS) $(IFLAGS) $(TESTFLAGS) $< -o $@.d
build/%.o: src/cpp/lib/%.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o $@ $<
	$(CXX) -MT $@ -MM $(CXXFLAGS) $(IFLAGS) $< -o $@.d
build/%.o: src/c/%.c 
	$(CC)  $(CFLAGS) $(IFLAGS) -c -o $@ $<
	$(CC) -MT $@ -MM $(CFLAGS) $(IFLAGS) $< -o $@.d

bin/gvcfgenotyper: src/cpp/gvcfgenotyper.cpp build/version.hh $(OBJS) $(HTSLIB)
	$(CXX) $(CXXFLAGS) -o $@   src/cpp/gvcfgenotyper.cpp $(OBJS) $(HTSLIB) $(IFLAGS) $(LFLAGS)
bin/test_gvcfgenotyper: build/version.hh $(OBJS) $(TESTOBJS) $(HTSLIB) build/gtest.a build/gtest_main.a
	$(CXX) $(CXXFLAGS) $(TESTFLAGS) -o $@ $(TESTOBJS) $(OBJS) $(IFLAGS) $(HTSLIB) $(LFLAGS) build/gtest.a build/gtest_main.a
.PHONY: test
test: bin/test_gvcfgenotyper bin/gvcfgenotyper
	bin/test_gvcfgenotyper
	bash -e src/bash/run_smoke_tests.sh
.PHONY: clean
clean:
	rm -rf build/* bin/*
