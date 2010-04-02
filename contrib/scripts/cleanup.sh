#! /bin/sh

if [ ! -d R ] ; then
  echo Only run cleanup from base dir
fi

cd contrib/bin
make clean
cd ../..

rm -rvf ..Rcheck/
rm tests/junk*
for x in cmake_install.cmake CTestTestfile.cmake CMakeCache.txt Makefile
do
  rm -vf $x
  find . -name $x -exec rm -v \{\} \;
done

find . -name CMakeFiles -exec rm -rvf \{\} \;
find . -name Testing -exec rm -rvf \{\} \;
find . -name *.so -exec rm -v \{\} \;
find . -name *.o -exec rm -v \{\} \;
find . -name *.a -exec rm -v \{\} \;
find . -name *dump -exec rm -v \{\} \;
find . -name *.aux -exec rm -v \{\} \;
find . -name *.dvi -exec rm -v \{\} \;
find . -name *.log -exec rm -v \{\} \;
find . -name *.orig -exec rm -v \{\} \;
find src/ -name *.dll -exec rm -v \{\} \;

