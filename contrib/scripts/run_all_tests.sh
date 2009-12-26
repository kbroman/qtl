#! /bin/sh

path=$1

echo "Run all functional/regression/unit/check tests for R/qtl on $path"
if [ ! -z $path -a -d $path ]; then
  cd $path
fi
echo -n "Using: "
pwd
if [ ! -d "contrib" ]; then
  echo "Incorrect path for R/qtl source"
  exit 1
fi
cwd=`pwd`

echo "Run the standard MQM regression tests"
cd $cwd
cd contrib/bin
rm CMakeCache.txt
cmake .
make
make testall
if [ "$?" -ne "0" ]; then
  echo "Test failed"
  exit 1
fi

echo "R CMD check"
cd $cwd
R CMD check .
if [ "$?" -ne "0" ]; then
  echo "Test failed"
  exit 1
fi

