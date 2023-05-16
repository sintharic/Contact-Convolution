SHELL = /usr/bin/env bash
default_target: all

COMPILER := g++
CPPFLAGS := -O2 -std=c++17
#WFLAGS   := -Wno-unused-command-line-argument -Wno-shift-count-overflow -Wno-unused-result
OBJECTS  := io.o sim.o
opts     :=
INC_PATH := /usr/local/include

io.o: header.h io.cpp
	@echo " compiling io.o ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) -c io.cpp -o io.o > error.log 2>&1

elast.o: header.h elast.cpp
	@echo " compiling elast.o ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) -c elast.cpp -o elast.o > error.log 2>&1

sim.o: header.h sim.cpp
	@echo " compiling sim.o ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) $(opts) -c sim.cpp -o sim.o > error.log 2>&1

all: $(OBJECTS)
	@echo " compiling sim.exe ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) $(OBJECTS) -o sim.exe > error.log 2>&1

sim: io.o
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) sim.cpp -o sim.exe > error.log 2>&1

GF:
	@echo " compiling GF.exe ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) GF.cpp -o GF.exe

cln:
	@rm -f $(OBJECTS)
	@rm -f sim.exe GF.exe