# Bayes++ GCC Makefile
# Bayesian filtering library

BOOST ?= ../../boost_1_29_0
INC = -I.. -I${BOOST}
CPPFLAGS = -O3 -DNDEBUG -Wall $(INC)

OBJS=bayesFlt.o bayesFltAlg.o matSup.o UdU.o covFlt.o infFlt.o infRtFlt.o itrFlt.o SIRFlt.o UDFlt.o unsFlt.o CIFlt.o

all : bayesFilter_gcc.lib

bayesFilter_gcc.lib : $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)


$(OBJS) : bayesFlt.h matSup.h uBLASmatrix.h
covFlt.o : covFlt.h
infFlt.o : infFlt.h
infRtFlt.o : infRtFlt.h uLAPACK.h
itrFlt.o : itrFlt.h
SIRFlt.o : SIRFlt.h
UDFlt.o : UDFlt.h
unsFlt.o : unsFlt.h
CIFlt.o : CIFlt.h
