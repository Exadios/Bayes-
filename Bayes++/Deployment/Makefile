# Bayes++ Project Managment Makefile
# Provides interface project development and sourceforge maintance
# Use bjam to compile and build Bayes++

RELEASE=2003.8-1

help:
	echo Usage: make [uploaddoc] [uploadstyle]

SFWEBLIST = *.htm Simple/simpleExample.cpp "ClassDocumentation"
uploaddoc:
	cd ..; rsync --verbose --recursive --relative --compress --times --rsh=ssh $(SFWEBLIST)  mistevens@shell.sf.net:/home/groups/b/ba/bayesclasses/htdocs


uploadstyle:
	rsync --verbose --relative --compress --times --rsh=ssh paper.css paper.png  mistevens@shell.sf.net:/home/groups/b/ba/bayesclasses/htdocs/style

archivedoc:
	cd ..; tar --create --file=bf-classdocumentation-$(RELEASE).tar.bz2 --use-compress-program=bzip2 ClassDocumentation/*
	cd ..; zip -r bf-classdocumentation-$(RELEASE).zip ClassDocumentation

archivesrc:
	cd ..; zip -r bf-C++source-$(RELEASE).zip BayesFilter Simple PV QuadCalib SLAM Test VClib  Jamfile Jamrules Jamfile.v2 projec-root.jam Makefile Bayes++.doxygen  Bayes++.sln .project Bayes++.htm "Bayesian Filtering Classes.htm"

buildallCVS :
	bjam -sBOOST=../boost/boost_CVS debug release compress_debug compress_release gcc-3.3.1 gcc-3.2 gcc-2.95 intellinux-70 msvcms . rtheta SLAM