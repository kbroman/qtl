echo "Starting regression tests"
cd test
dir
PATH=%PATH%;c:\Program Files\R\R-2.9.1\bin
# valgrind ../src/mqm -T=0 -V > ../test/MQM_test0.txt
..\mqm -T=0 -V > MQM_test0.txt

echo "Finalized regression tests"

