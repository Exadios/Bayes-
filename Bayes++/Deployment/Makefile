# Bayes++ Project Managment Makefile
#   Provides interface for project development and sourceforge maintainance
#   Use bjam to compile and build Bayes++

RELEASE=2005
RELEASE_SUB=-1
BOOST_VERSION=boost_1_32_0

help:
	echo Usage: make [uploaddoc] [uploadstyle] [archivedoc] [archivesrc] [buildallCVS]

SFWEBLIST = *.html Simple/simpleExample.cpp "ClassDocumentation"
uploaddoc:
	cd ..; rsync --verbose --recursive --relative --compress --times --rsh=ssh $(SFWEBLIST)  mistevens@shell.sf.net:/home/groups/b/ba/bayesclasses/htdocs

uploadstyle:
	rsync --verbose --relative --compress --times --rsh=ssh paper.css paper.png  mistevens@shell.sf.net:/home/groups/b/ba/bayesclasses/htdocs/Deployment


archivedoc:
	cd ..; tar --create --file=Deployment/bf-classdocumentation-$(RELEASE).tar.bz2 --use-compress-program=bzip2 ClassDocumentation/*
	cd ..; zip -r Deployment/bf-classdocumentation-$(RELEASE).zip ClassDocumentation

archivesrc:
	ln -s "../../Bayes++$(RELEASE)" Bayes++
	zip -r "Bayes++/bf-C++source-$(RELEASE)$(RELEASE_SUB).zip" Bayes++/BayesFilter Bayes++/Simple Bayes++/PV Bayes++/QuadCalib Bayes++/SLAM Bayes++/Test Bayes++/VClib  Bayes++/Jamfile.v2 Bayes++/project-root.jam Bayes++/Makefile Bayes++/Bayes++.sln Bayes++/.project Bayes++/Bayes++.doxygen Bayes++/Bayes++.html "Bayes++/Bayesian Filtering Classes.html"  -x '*/CVS/*'
	rm Bayes++

buildall :
	cd ..; bjam -sBOOST_ROOT=../Boost/$(BOOST_VERSION) debug release compress_debug compress_release gcc-3.3 gcc-3.4 intel-linux-8.0 msvc-7.1/runtime-link=static  . SLAM rtheta
	cd ..; bjam -sBOOST_ROOT=../Boost/$(BOOST_VERSION) debug release msvc-7.0 . SLAM rtheta

buildallCVS :
	cd ..; bjam -sBOOST_ROOT=../Boost/boost_CVS debug release compress_debug compress_release gcc-3.3 gcc-3.4 intel-linux-8.0 msvc-7.1 . SLAM rtheta
	cd ..; bjam -sBOOST_ROOT=../Boost/boost_CVS debug release msvc-7.0 . SLAM rtheta
	