#! /bin/sh
#
# Usage: contrib/scripts/check_rqtl.sh

path=$1

echo "* Run all R CMD check for R/qtl on $path"
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


echo "* Run R CMD check"
cd $cwd
sh contrib/scripts/cleanup.sh
R CMD check .
if [ "$?" -ne "0" ]; then
  echo "Test 'R CMD check' failed"
  exit 1
fi


