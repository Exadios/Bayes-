# Bayes++ Project Makefile
#
# The BayesFilter library and examples as built using the Boost Build system.
# The Boost build system provides a flexible and toolset independant replacment for make.
# Bayes++ also uses several Boost template libraries. Therefore it is simple
# to combine the build process for Bayes++, with the libraries it requires.
#
# The Boost Build system uses the 'bjam' tool.
# This Makefile just wraps up the bjam command for those used to make.
# Use make ARGS='aa'  to pass additional arguments to bjam

# Modify this variable to locate bjam
BJAM ?= bjam

GOBJAM = $(BJAM) $(ARGS)

# Build everything
buildAll:
	$(GOBJAM)


# Build individual examples
buildSimple:
	$(GOBJAM) '<Simple>simple'

buildPV:
	$(GOBJAM) '<PV>pv'

buildQuadCalib:
	$(GOBJAM) '<QuadCalib>quadcalib'

