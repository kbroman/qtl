#! /bin/bash

echo "Starting profiling"
cd test
gprof ../src/sMQM -- -T=0 -V
echo "Finalized profiling"
