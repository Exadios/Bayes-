# Bayes++ Project Makefile
# 
# Bayes++ uses several Boost template libraries.
# The Boost.build system (and bjam) is also used to create the Bayes++ examples.
# It provides a flexible and toolset independant system to build the examples and
# a BayesFilter static library.
#
# This Makefile just wraps up the bjam command for those used to make.

# Modify these settings to locate Boost and bjam
BOOST_ROOT ?= ../boost_1_29_0
BJAM ?= ${BOOST_ROOT}/bjam
TOOL ?= gcc
BJAMFLAGS ?= -sALL_LOCATE_TARGET=bin

# Build everything
bayes++:
	${BJAM} -sBOOST_ROOT=$(BOOST_ROOT) ${BJAMFLAGS}


# Build individual examples
simple:
	${BJAM} -sBOOST_ROOT=$(BOOST_ROOT) ${BJAMFLAGS} '<Simple>simple'

rthetatest:
	${BJAM} -sBOOST_ROOT=$(BOOST_ROOT) ${BJAMFLAGS} '<rtheta>rtheta'
