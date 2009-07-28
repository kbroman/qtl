cmake_minimum_required(VERSION 2.6)

IF (NOT CMAKE_MODULE_PATH)
  set(CMAKE_MODULE_PATH ".")
ENDIF()

SET(CMAKE_BUILD_TYPE Debug)

IF ( CMAKE_BUILD_TYPE MATCHES Debug )
  ADD_DEFINITIONS(-Wall -D_DEBUG)
ENDIF ( CMAKE_BUILD_TYPE MATCHES Debug )

ENABLE_TESTING()

FIND_PACKAGE(RLibs)

include_directories (../../src)
add_definitions(-DSTANDALONE)
add_definitions(-pg)

set(LINK_FLAGS "-pg")

link_libraries(${R_LIBRARY})

add_executable (mqm 
../../src/mqmaugment.cpp    
../../src/mqmeliminate.h  
../../src/mqmmapqtl.h     
../../src/mqmregression.cpp
../../src/mqmaugment.h      
../../src/mqmextra.cpp    
../../src/mqmmixture.cpp  
../../src/mqmregression.h
../../src/mqmdatatypes.cpp  
../../src/mqmextra.h      
../../src/mqmmixture.h    
../../src/mqmscan.cpp
../../src/mqmdatatypes.h    
../../src/mqm.h           
../../src/mqmprob.cpp     
../../src/mqmscan.h
../../src/mqmeliminate.cpp  
../../src/mqmmapqtl.cpp   
../../src/mqmprob.h       
../../src/standalone.h
../../src/util.h
../../src/util.c
mqmmain.cpp
)

ADD_TEST(TestStandaloneMQM ./scripts/regression_tests.sh)

