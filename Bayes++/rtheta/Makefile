# Bayes++ GCC Makefile
# Range Angle test for all filter schemse

UBLASLIB ?= ..
BOOSTLIB ?= ../../boost_1_28_0
INC = -I.. -I../uBLASmatSup -I${UBLASLIB} -I${BOOSTLIB}
CC=g++
CPPFLAGS = -O3 -DNDEBUG -finline-limit=1000 $(INC)

OBJS=rtheta.o test.o

all : rtheta

rtheta : $(OBJS) ../BayesFilter/bayesFilter_gcc.lib

