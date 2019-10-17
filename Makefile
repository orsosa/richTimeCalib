.DELETE_ON_ERROR:
.PHONY: all

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

ifndef HIPO4INC
$(error HIPO4INC is not set)
endif

ifndef HIPO4LIB
$(error HIPO4LIB is not set)
endif

ifeq ($(wildcard .obj/.),)
	@echo ".obj not found, creating the directory..."
	mkdir .obj
endif
ifeq ($(wildcard bin/.),)
	@echo "bin not found, creating the directory..."
	mkdir bin
endif

HIPOCFLAGS  := -I$(HIPO4INC)
HIPOLIBS    := -L$(HIPO4LIB) -lhipo4

LZ4LIBS     := -L$(LZ4DIR)/lib -llz4
LZ4INCLUDES := -I$(LZ4DIR)/lib


INC_DIR   := -Iinclude 

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS) -std=c++11
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

AR	  = ar
ARFLAGS	  = -cvr #create,verbose,quick (don't check for replacement, otherwise use r instead)


all:  bin/richTiming bin/richTimeOffsets bin/richTimeWalks bin/richTimeCheck


bin/%: .obj/%.o
	$(CXX) -o $@ $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS)  $(ROOTLIBS) 

.obj/%.o: src/%.cxx
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES) $(INC_DIR)

install:
	chmod +x bin/*
	chmod +x script/*

clean:
	rm bin/* 
