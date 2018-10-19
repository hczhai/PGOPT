SHELL = /bin/bash
OST := $(shell uname -s)
HOME := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

CC = gcc
FPC = f2py
FPCFLAGS = --fcompiler=intelem --f77exec=ifort --f90exec=ifort -lmkl_core -lmkl_intel_lp64 -lmkl_sequential

SRCDIR   := $(HOME)
SRCFILES = $(shell find $(SRCDIR) -name "*.f90")
OBJFILES = $(addsuffix .so, $(basename $(notdir $(SRCFILES))))

.PHONY : all clean
default_target: all

all: $(OBJFILES)
	@echo 'compile completed.'

clean:
	rm -rf $(HOME)/*.so
	rm -rf $(HOME)/*.so.dSYM
	rm -rf $(HOME)/intel

%.so: %.f90
ifeq ($(OST),Darwin)
	$(FPC) -c -m $(basename $(notdir $<)) $< $(FPCFLAGS)
	install_name_tool -change libmkl_core.dylib $(MKLROOT)/lib/libmkl_core.dylib $@
	install_name_tool -change libmkl_intel_lp64.dylib $(MKLROOT)/lib/libmkl_intel_lp64.dylib $@
	install_name_tool -change libmkl_sequential.dylib $(MKLROOT)/lib/libmkl_sequential.dylib $@
else
	$(FPC) -c -m $(basename $(notdir $<)) $< $(FPCFLAGS)
endif
