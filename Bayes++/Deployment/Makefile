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
	cd ../..; zip -r Bayes++/bf-C++source-$(RELEASE).zip Bayes++/BayesFilter Bayes++/Simple Bayes++/PV Bayes++/QuadCalib Bayes++/SLAM Bayes++/Test Bayes++/VClib  Bayes++/Jamfile Bayes++/Jamrules Bayes++/Jamfile.v2 Bayes++/project-root.jam Bayes++/Makefile Bayes++/Bayes++.sln Bayes++/.project Bayes++/Bayes++.doxygen Bayes++/Bayes++.htm "Bayes++/Bayesian Filtering Classes.htm"

buildallCVS :
	cd ..; bjam -sBOOST_ROOT=../boost/boost_CVS debug release compress_debug compress_release gcc-3.3 gccms-3.4 gccms-2.95 intellinux-70 msvcms . SLAM rtheta
