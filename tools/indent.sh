#!/bin/bash
#
# Install emacs: 
#	sudo apt-get install emacs
#
# Get the indentation emacs script as follows: 
#  	wget https://raw.github.com/Glavin001/atom-beautify/master/src/beautifiers/fortran-beautifier/emacs-fortran-formating-script.lisp#L44
#
# To use it just call this script from the main LATTE folder. 
# 	./tools/indent.sh

LATTE_PATH=$HOME/LATTE/

for file in $LATTE_PATH/src/*.f90
do 
  emacs -batch -l $LATTE_PATH/tools/emacs-fortran-formating-script.lisp -f f90-batch-indent-region $file
done	


