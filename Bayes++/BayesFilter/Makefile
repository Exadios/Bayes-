# Bayes++ GCC Makefile
# Bayesian filtering library

UBLASLIB ?= ..
BOOSTLIB ?= ../../boost_1_28_0
INC = -I.. -I../uBLASmatSup -I${UBLASLIB} -I${BOOSTLIB}
CPPFLAGS = -O3 -DNDEBUG -Wall $(INC)

OBJS=bayesFlt.o bayesFltAlg.o matSup.o UdU.o covFlt.o infFlt.o infRtFlt.o itrFlt.o SIRFlt.o UDFlt.o unsFlt.o

all : bayesFilter_gcc.lib

bayesFilter_gcc.lib : $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)


$(OBJS) : bayesFlt.h matSup.h ../uBLASmatSup/uBLASmatrix.h
covFlt.o : covFlt.h
infFlt.o : infFlt.h
infRtFlt.o : infRtFlt.h ../uBLASmatSup/uLAPACK.h
itrFlt.o : itrFlt.h
SIRFlt.o : SIRFlt.h
UDFlt.o : UDFlt.h
unsFlt.o : unsFlt.h
