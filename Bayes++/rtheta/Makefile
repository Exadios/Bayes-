# Bayes++ rtheta Makefile
# Range Angle test for all filter schemes
#

# Modify these settings to locate bjam
BJAM = bjam

GOBJAM = $(BJAM) $(ARGS)

# Build everything
buildRtheta:
	$(GOBJAM)

buildRtheta-boost-CVS-debug:
	$(GOBJAM) -s"BUILD=boost-CVS-debug"

buildRtheta-boost-CVS-release:
	$(GOBJAM) -s"BUILD=boost-CVS-release"
