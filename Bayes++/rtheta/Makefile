# Bayes++ rtheta Makefile
# Range Angle test for all filter schemse
#

# Modify these settings to locate Boost and bjam
BOOST_ROOT = ../../boost
BJAM = bjam
ALL_LOCATE_TARGET = bin

GOBJAM = $(BJAM) $(ARGS) -sBOOST_ROOT=$(BOOST_ROOT) -sALL_LOCATE_TARGET=$(ALL_LOCATE_TARGET)

# Build everything
buildRtheta:
	$(GOBJAM)
