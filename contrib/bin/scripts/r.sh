#!/bin/sh

echo Testing $1/$2
cd $1
R --no-save --no-restore --no-readline --slave < $2
