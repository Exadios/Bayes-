# Bayes++ GCC Makefile
# Simple filter example

BOOST ?= ../../boost_1_29_0
INC = -I.. -I../uBLASmatSup -I${BOOST}
CPPFLAGS = -DNDEBUG -Wall $(INC)

all : simpleFilter

simpleFilter : simpleFilter.cpp ../BayesFilter/bayesFilter_gcc.lib

