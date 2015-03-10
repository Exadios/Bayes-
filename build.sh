#!/bin/bash

bjam toolset=gcc link=static cxxflags=-Wno-unused-local-typedefs address-model=64 -sBOOST_ROOT="../boost_1_57_0"
