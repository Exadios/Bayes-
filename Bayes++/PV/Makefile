# Bayes++ GCC Makefile
# PV - Position Velocity filter

BOOST ?= ../../boost_1_29_0
INC = -I.. -I../uBLASmatSup -I${BOOST}
CPPFLAGS = -DNDEBUG -Wall $(INC)

all : PV

PV : PV.cpp ../BayesFilter/bayesFilter_gcc.lib

