# Bayes++ GCC3 Makefile
# PV - Position Velocity filter

UBLASLIB ?= ..
BOOSTLIB ?= ../../boost_1_28_0
INC=-I.. -I../UBLASmatSup -I${UBLASLIB} -I${BOOSTLIB}
CC=g++
CPPFLAGS=-DNDEBUG -Wall $(INC)

OBJS=PV.o

all : PV

PV : $(OBJS) ../BayesFilter/bayesFilter_gcc3.lib

