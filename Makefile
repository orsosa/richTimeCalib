.DELETE_ON_ERROR:
.PHONY: all checkdir

ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

ifndef CLAS12TOOL
$(error CLAS12TOOL is not set)
endif

HIPO4INC := $(if $(HIPO4INC),$(HIPO4INC),$(CLAS12TOOL)/Hipo4 )
HIPO4LIB := $(if $(HIPO4LIB),$(HIPO4LIB),$(CLAS12TOOL)/lib )

HIPOCFLAGS  := -I$(HIPO4INC)
HIPOLIBS    := -L$(HIPO4LIB) -lhipo4

LZ4LIBS     := -L$(LZ4DIR)/lib -llz4
LZ4INCLUDES := -I$(LZ4DIR)/lib

CLAS12ANALIB  := -L$(CLAS12ANAHIPO4)/shlib 
INC_DIR   := -I$(CLAS12ANAHIPO4)/include -I$(HIPO4ROOT)/Schemas -Iinclude

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS) -std=c++11
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

AR	  = ar
ARFLAGS	  = -cvr #create,verbose,quick (don't check for replacement, otherwise use r instead)


all:  checkdir bin/richTiming_app bin/richTimeOffsets bin/richTimeWalks bin/richTimeCheck

#richTiming:  richTiming

#richTiming_expert:  richTiming_expert

checkdir:
ifeq ($(wildcard .obj/.),)
	@echo ".obj not found, creating the directory..."
	mkdir .obj
endif
ifeq ($(wildcard bin/.),)
	@echo "bin not found, creating the directory..."
	mkdir bin
endif


bin/%: .obj/%.o
	$(CXX) -o $@ $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS)  $(ROOTLIBS) $(CLAS12ANALIB) 


.obj/%.o: src/%.cxx
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES) $(INC_DIR)

clean:
	rm richTiming richTiming_expert richTiming_app
