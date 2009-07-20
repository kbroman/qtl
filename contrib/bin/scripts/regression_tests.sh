#! /bin/bash

echo "Starting regression tests"
cd test
ls -l
# valgrind ../src/mqm -T=0 -V > ../test/MQM_test0.txt
../mqm -T=0 -V > MQM_test0.txt

echo "Finalized regression tests"
exit 0
