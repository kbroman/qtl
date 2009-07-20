cmake_minimum_required(VERSION 2.6)

IF (NOT CMAKE_MODULE_PATH)
  set(CMAKE_MODULE_PATH "../tools/cmake-support/modules")
ENDIF()

SET(CMAKE_BUILD_TYPE Debug)
ENABLE_TESTING()

FIND_PACKAGE(RLibs)

include_directories (.)
add_definitions(-DSTANDALONE)
add_definitions(-pg)

set(LINK_FLAGS "-pg")

link_libraries(${R_LIBRARY})

add_executable (mqm 
mqmaugment.cpp    mqmeliminate.h  mqmmapqtl.h     mqmregression.cpp
mqmaugment.h      mqmextra.cpp    mqmmixture.cpp  mqmregression.h
mqmdatatypes.cpp  mqmextra.h      mqmmixture.h    mqmscan.cpp
mqmdatatypes.h    mqm.h           mqmprob.cpp     mqmscan.h
mqmeliminate.cpp  mqmmapqtl.cpp   mqmprob.h       mqmmain.cpp
standalone.h
util.h
util.c
)

ADD_TEST(TestStandaloneMQM ../scripts/regression_tests.sh)

