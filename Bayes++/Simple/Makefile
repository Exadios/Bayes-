# Bayes++ GCC3 Makefile
# Simple Example

UBLASLIB ?= ..
BOOSTLIB ?= ../../boost_1_28_0
INC=-I.. -I../UBLASmatSup -I${UBLASLIB} -I${BOOSTLIB}
CC=g++
CPPFLAGS=-DNDEBUG -Wall $(INC)

OBJS=simpleExample.o

all : simpleExample

simpleExample : $(OBJS) ../BayesFilter/bayesFilter_gcc3.lib

