# Bayes++ GCC Makefile
# Range Angle test for all filter schemse

TESTING ?= ../../libs
BOOST ?= ../../boost_1_29_0
INC=-I.. -I../uBLASmatSup -I${BOOST} -I$(TESTING)
CC=g77	# g77 Required to link correctly with lapack which requires fortran library but not startup code
CPPFLAGS=-O3 -DNDEBUG $(INC)
LDFLAGS=-llapack -lstdc++

OBJS=rtheta.o tests.o

all : rtheta

rtheta : $(OBJS) ../BayesFilter/bayesFilter_gcc.lib

