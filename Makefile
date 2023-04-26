SHELL = /usr/bin/env bash
default_target: all

COMPILER := g++
CPPFLAGS := -O2 -std=c++11
#WFLAGS   := -Wno-unused-command-line-argument -Wno-shift-count-overflow -Wno-unused-result
OBJECTS  := io.o sim.o
opts     :=

io.o: header.h io.cpp
	@echo " compiling io.o ..."
	@$(COMPILER) $(CPPFLAGS) -c io.cpp -o io.o > error.log 2>&1

sim.o: header.h sim.cpp
	@echo " compiling sim.o ..."
	@$(COMPILER) $(CPPFLAGS) $(opts) -c sim.cpp -o sim.o > error.log 2>&1

all: $(OBJECTS)
	@echo " compiling sim.exe ..."
	@$(COMPILER) $(CPPFLAGS) $(OBJECTS) -o sim.exe > error.log 2>&1

GF:
	@echo " compiling GF.exe ..."
	@$(COMPILER) $(CPPFLAGS) GF.cpp -o GF.exe

cln:
	@rm -f $(OBJECTS)
	@rm -f sim.exe GF.exe