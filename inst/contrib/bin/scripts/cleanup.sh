#! /bin/bash

for x in cmake_install.cmake CTestTestfile.cmake CMakeCache.txt Makefile
do
  rm -vf $x
  find . -name $x -exec rm -v \{\} \;
done

find . -name CMakeFiles -exec rm -rvf \{\} \;
find . -name Testing -exec rm -rvf \{\} \;
find . -name *.so -exec rm -v \{\} \;
find . -name *.dll -exec rm -v \{\} \;
find . -name *.o -exec rm -v \{\} \;
find . -name *.a -exec rm -v \{\} \;
find . -name *dump -exec rm -v \{\} \;

