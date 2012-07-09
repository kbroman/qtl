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
find . -name *.eps -exec rm -v \{\} \;
find . -name *~ -exec rm -v \{\} \;
find . -name Rplots.pdf -exec rm -v \{\} \;
find . -name *.Rout -exec rm -v \{\} \;
find src/ -name *.dll -exec rm -v \{\} \;

rm inst/docs/Sources/MQM/MQM-tour.tex
mv inst/docs/Sources/MQM/MQM-tour.R inst/docs/
mv inst/docs/Sources/MQM/MQM-tour.pdf inst/docs/

rm -rvf src-i386/
rm -rvf src-x86_64/
