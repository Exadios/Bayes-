# Bayes++ Project Managment Makefile
# Provides interface project development and sourceforge maintance
# Use bjam to compile and build Bayes++

help:
	echo Usage: make uploadweb

SFWEBLIST = *.htm Simple/simpleExample.cpp "ClassDocumentation"
uploadweb:
	rsync --verbose --recursive --relative --compress --times --rsh=ssh $(SFWEBLIST)  mistevens@shell.sf.net:/home/groups/b/ba/bayesclasses/htdocs
