#! /bin/bash
#
# Create a diff from another repostitory (this one should be standalone 
# and the other a master).

if [ ! -d .git ]; then
  echo Should be in root of repo
  exit 1
fi

git checkout standalone
cd src
ls --color=never *.cpp *.h > standalone.lst
git checkout master
cat standalone.lst | grep -v main | xargs git diff standalone > standalone.patch

