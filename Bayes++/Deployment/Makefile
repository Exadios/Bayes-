# Bayes++ Project Managment Makefile
# Provides interface project development and sourceforge maintance
# Use bjam to compile and build Bayes++

RELEASE=2004.2.1
BOOST_VERSION=boost_1_31_0

help:
	echo Usage: make [uploaddoc] [uploadstyle] [archivedoc] [archivesrc] [buildallCVS]

SFWEBLIST = *.html Simple/simpleExample.cpp "ClassDocumentation"
uploaddoc:
	cd ..; rsync --verbose --recursive --relative --compress --times --rsh=ssh $(SFWEBLIST)  mistevens@shell.sf.net:/home/groups/b/ba/bayesclasses/htdocs

uploadstyle:
	rsync --verbose --relative --compress --times --rsh=ssh paper.css paper.png  mistevens@shell.sf.net:/home/groups/b/ba/bayesclasses/htdocs/Deployment


archivedoc:
	cd ..; tar --create --file=bf-classdocumentation-$(RELEASE).tar.bz2 --use-compress-program=bzip2 ClassDocumentation/*
	cd ..; zip -r bf-classdocumentation-$(RELEASE).zip ClassDocumentation

archivesrc:
	cd ../..; zip -r Bayes++/bf-C++source-$(RELEASE).zip Bayes++/BayesFilter Bayes++/Simple Bayes++/PV Bayes++/QuadCalib Bayes++/SLAM Bayes++/Test Bayes++/VClib  Bayes++/Jamfile Bayes++/Jamrules Bayes++/Jamfile.v2 Bayes++/project-root.jam Bayes++/Makefile Bayes++/Bayes++.sln Bayes++/.project Bayes++/Bayes++.doxygen Bayes++/Bayes++.htm "Bayes++/Bayesian Filtering Classes.htm"

buildall :
	cd ..; bjam -sBOOST_ROOT=../boost/$(BOOST_VERSION) debug release compress_debug compress_release gcc intellinux-70 msvc-7 . SLAM rtheta
	cd ..; bjam -sBOOST_ROOT=../boost/$(BOOST_VERSION) debug release gcc-2.95 . SLAM

buildallCVS :
	cd ..; bjam -sBOOST_ROOT=../boost/boost_CVS debug release compress_debug compress_release gcc intellinux-70 msvc-7 . SLAM rtheta
	cd ..; bjam -sBOOST_ROOT=../boost/boost_CVS debug release gcc-2.95 . SLAM
	