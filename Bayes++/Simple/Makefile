# Bayes++ GCC Makefile
# PV - Position Velocity filter

UBLASLIB ?= ..
BOOSTLIB ?= ../../boost_1_28_0
INC = -I.. -I../uBLASmatSup -I${UBLASLIB} -I${BOOSTLIB}
CPPFLAGS = -DNDEBUG -Wall $(INC)

all : simpleFilter

simpleFilter : simpleFilter.cpp ../BayesFilter/bayesFilter_gcc.lib

