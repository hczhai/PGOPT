SHELL = /bin/bash
OST := $(shell uname -s)
HOME := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))

CC = gcc
FPC = f2py
FPCFLAGS = --fcompiler=gfortran

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
	$(FPC) -c -m $(basename $(notdir $<)) $< $(FPCFLAGS)
