# Bayes++ Project Makefile
# 
# Bayes++ uses several Boost template libraries.
# The Boost.build system (and bjam) is also used to create the Bayes++ examples.
# It provides a flexible and toolset independant system to build the examples and
# a BayesFilter static library.
#
# This Makefile just wraps up the bjam command for those used to make.
# Use ARGS='xx' to pass arguments to bjam

# Modify these settings to locate Boost and bjam
BOOST_ROOT ?= ../boost_1_29_0
BJAM ?= ${BOOST_ROOT}/bjam
TOOL ?= gcc
ALL_LOCATE_TARGET ?= bin

GOBJAM = ${BJAM} -sBOOST_ROOT=$(BOOST_ROOT) -sALL_LOCATE_TARGET=$(ALL_LOCATE_TARGET)  $(ARGS)

# Build everything
all:
	${GOBJAM}


# Build individual examples
Simpletest:
	${GOBJAM} '<Simple>simple'

PVtest:
	${GOBJAM} '<PV>PV'

QuadCalibtest:
	${GOBJAM} '<QuadCalib>Quadcalib'

RThetatest:
	${GOBJAM} '<rtheta>rtheta'

FastSLAMtest:
	${GOBJAM} '<SLAM>testfastslam'
