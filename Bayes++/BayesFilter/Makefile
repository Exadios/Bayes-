# Bayes++ GCC3 Makefile
# Bayesian filtering library

UBLASLIB ?= ..
BOOSTLIB ?= ../../boost_1_28_0
INC=-I.. -I../UBLASmatSup -I${UBLASLIB} -I${BOOSTLIB}
CPPFLAGS=-O3 -DNDEBUG -finline-limit=1000 -Wall $(INC)

OBJS=bayesFlt.o bayesFltAlg.o covFlt.o infFlt.o infRtFlt.o itrFlt.o matSup.o SIRFlt.o UDFlt.o UdU.o unsFlt.o

all : bayesFilter_gcc3.lib

bayesFilter_gcc3.lib : $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)


$(OBJS) : bayesFlt.h matSup.h
covFlt.o : covFlt.h
infFlt.o : infFlt.h
infRtFlt.o : infRtFlt.h
itrFlt.o : itrFlt.h

