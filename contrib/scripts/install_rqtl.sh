#! /bin/sh
#
# Usage: contrib/scripts/install_rqtl.

path=$1

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

sh contrib/scripts/cleanup.sh

rqtl_version=`grep Version DESCRIPTION | awk '{ print $2; }'`

echo "* Run R CMD check"
cd $cwd
cd ..
R CMD build $cwd
R CMD INSTALL qtl_${rqtl_version}.tar.gz
rm -v qtl_${rqtl_version}

cd $cwd
time R --no-save --no-restore --no-readline --slave < ./tests/test_qtl.R 

