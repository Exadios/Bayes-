# Bayes++ GCC3 Makefile
# FastSLAM Test

UBLASLIB ?= ../../ublas
BOOSTLIB ?= ../../boost_1_28_0
INC=-I.. -I../UBLASmatSup -I${UBLASLIB} -I${BOOSTLIB}
CC=g++
CPPFLAGS=-O3 -DNDEBUG -finline-limit=1000 $(INC)

OBJS=testFastSLAM.o fastSLAM.o kalmanSLAM.o

all : testFastSLAM 

testFastSLAM : $(OBJS) ../BayesFilter/bayesFilter_gcc3.lib

