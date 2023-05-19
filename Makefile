SHELL = /usr/bin/env bash
default_target: test

COMPILER := g++
CPPFLAGS := -O2 -std=c++17
#WFLAGS   := -Wno-unused-command-line-argument -Wno-shift-count-overflow -Wno-unused-result
OBJECTS  := io.o ElasticBody.o Indenter.o Interaction.o
opts     :=
INC_PATH := /usr/local/include

io.o: header.h io.cpp
	@echo " compiling io.o ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) -c io.cpp -o io.o > error.log 2>&1

ElasticBody.o: header.h ElasticBody.h ElasticBody.cpp
	@echo " compiling ElasticBody.o ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) -c ElasticBody.cpp -o ElasticBody.o > error.log 2>&1

Indenter.o: header.h Indenter.h Indenter.cpp
	@echo " compiling Indenter.o ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) -c Indenter.cpp -o Indenter.o > error.log 2>&1

Interaction.o: header.h ElasticBody.h Indenter.h Interaction.h Interaction.cpp
	@echo " compiling Interaction.o ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) -c Interaction.cpp -o Interaction.o > error.log 2>&1

sim.o: header.h ElasticBody.h Indenter.h Interaction.h sim.cpp
	@echo " compiling sim.o ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) $(opts) -c sim.cpp -o sim.o > error.log 2>&1

test.o: header.h ElasticBody.h Indenter.h Interaction.h test.cpp
	@echo " compiling test.o ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) $(opts) -c test.cpp -o test.o > error.log 2>&1

sim: $(OBJECTS) sim.o
	@echo " compiling sim.exe ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) $(OBJECTS) sim.o -o sim.exe > error.log 2>&1

test: $(OBJECTS) test.o
	@echo " compiling test.exe ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) $(OBJECTS) test.o -o test.exe > error.log 2>&1

GF:
	@echo " compiling GF.exe ..."
	@$(COMPILER) $(CPPFLAGS) -I $(INC_PATH) GF.cpp -o GF.exe

cln:
	@rm -f $(OBJECTS) sim.o test.o
	@rm -f sim.exe test.exe GF.exe