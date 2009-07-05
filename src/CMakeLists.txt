cmake_minimum_required(VERSION 2.6)

IF (NOT CMAKE_MODULE_PATH)
  set(CMAKE_MODULE_PATH "../cmake-modules/")
ENDIF()

ENABLE_TESTING()

FIND_PACKAGE(RLibs)

include_directories (.)
add_definitions(-DSTANDALONE)
add_definitions(-pg)

set(LINK_FLAGS "-pg")

link_libraries(${R_LIBRARY})

add_executable (sMQM 
MQM.h
MQMaugment.cpp
MQMaugment.h
MQMdatatypes.cpp
MQMdatatypes.h
MQMeliminate.cpp
MQMeliminate.h
MQMextra.cpp
MQMextra.h
MQMmain.cpp
MQMmapqtl.cpp
MQMmapqtl.h
MQMmixture.cpp
MQMmixture.h
MQMprob.cpp
MQMprob.h
MQMregression.cpp
MQMregression.h
MQMscan.cpp
MQMscan.h
standalone.h
util.h
util.c
)

ADD_TEST(TestStandaloneMQM ../scripts/regression_tests.sh)

