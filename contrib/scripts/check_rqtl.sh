#! /bin/sh

echo "* Running R CMD check $1"
cd $1
cd ..
R CMD check $1

