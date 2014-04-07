#!/bin/sh

export DOXYGEN_DIR="../doc/doxygen"

cd $DOXYGEN_DIR

# 1. make the PDF
cd latex
make all 

# 2. create shortcut to the PDF and the main HTML page
cd ..
/bin/ln -fs latex/refman.pdf ODTmanual.pdf
/bin/ln -fs html/index.html ODTmanual.html

