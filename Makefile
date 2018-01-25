.PHONY: all
all: bin/gvcfgenotyper bin/test_gvcfgenotyper

# hard-coded version
VERSION_MAJOR=19-01-2018
# if we are in a git repo and if git binary is available, add git hash + branch info
GIT_HASH = $(shell git rev-parse --abbrev-ref HEAD 2> /dev/null)_$(shell git describe --always 2> /dev/null)
ifneq "$(GIT_HASH)" "_"
	VERSION= -DGG_VERSION=\"$(VERSION_MAJOR)_$(GIT_HASH)\"
else
	VERSION= -DGG_VERSION=\"$(VERSION_MAJOR)\"
endif

CC=gcc
CXX=g++

CXXFLAGS= -std=c++11 -O2 $(VERSION)
CFLAGS = -O2 $(VERSION)


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

# spdlog
SPDLOGDIR=external/spdlog
IFLAGS += -I$(SPDLOGDIR)

#special builds for debugging and profiling
debug: CXXFLAGS = -std=c++11 -g -O1 -Wall $(VERSION)
debug: CFLAGS = -g -O1 $(VERSION)
debug: all

profile: CXXFLAGS = -std=c++11 -pg $(VERSION)
profile: CFLAGS =  -pg $(VERSION)
profile: all


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

bin/gvcfgenotyper: src/cpp/gvcfgenotyper.cpp  $(OBJS) $(HTSLIB)
	$(CXX) $(CXXFLAGS) -o $@   src/cpp/gvcfgenotyper.cpp $(OBJS) $(HTSLIB) $(IFLAGS) $(LFLAGS)
bin/test_gvcfgenotyper:  $(OBJS) $(TESTOBJS) $(HTSLIB) build/gtest.a build/gtest_main.a
	$(CXX) $(CXXFLAGS) $(TESTFLAGS) -o $@ $(TESTOBJS) $(OBJS) $(IFLAGS) $(HTSLIB) $(LFLAGS) build/gtest.a build/gtest_main.a
.PHONY: test
test: bin/test_gvcfgenotyper bin/gvcfgenotyper
	bin/test_gvcfgenotyper
	bash -e src/bash/run_smoke_tests.sh
.PHONY: clean
clean:
	rm -rf build/* bin/*
