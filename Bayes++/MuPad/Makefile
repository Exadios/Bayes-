# Make Bfilter for MuPAD Module Generator

INC=-I.. -I../UBLASmatSup
MODULES=$(PROGRAMFILES)/SciFace/MuPAD Pro 2.0/modules
MMG="$(MODULES)/bin/Mmg.exe"

# Run Module Generator
# options:
#  -b		batch mode
#  -a		modattr for unloading
#  -oc  	CL options:  use bigger compiler heap, enable exception handling and RTTI
OPTS=-b -a unload -oc "-Zm150 -GX -GR $(INC)"
DOPTS=-g -D_DEBUG -oc "-MDd"

all:
	echo "No Target (debug/release)"

debug : bfilterD.mdm
	copy bfilterD.mdm "$(MODULES)/bfilter.mdm"

release : bfilterR.mdm
	copy bfilterR.mdm "$(MODULES)/bfilter.mdm"

clean:
	del *.mdm

# Use release clapack as can't find debug symbols
bfilterD.mdm : bfilter.cpp MuPadConvert.h "../VClib/BayesFilter_uBLAS_MD7d.lib" Makefile
	$(MMG) $(OPTS) $(DOPTS) "../VClib/BayesFilter_uBLAS_MD7d.lib" "../../CLAPACK3/Debug/clapack7d.lib" bfilter.cpp
	move /Y bfilter.mdm bfilterD.mdm

bfilterR.mdm : bfilter.cpp MuPadConvert.h "../VClib/BayesFilter_uBLAS_MD7.lib" Makefile
	$(MMG) $(OPTS) "../VClib/BayesFilter_uBLAS_MD7.lib" "../../CLAPACK3/Release/clapack7.lib" bfilter.cpp
	move /Y bfilter.mdm bfilterR.mdm
