cmake_minimum_required(VERSION 2.6)

IF (NOT CMAKE_MODULE_PATH)
  set(CMAKE_MODULE_PATH ".")
ENDIF()

SET(CMAKE_BUILD_TYPE Debug)

IF ( CMAKE_BUILD_TYPE MATCHES Debug )
  ADD_DEFINITIONS(-Wall -D_DEBUG)
ENDIF ( CMAKE_BUILD_TYPE MATCHES Debug )

ENABLE_TESTING()
ADD_CUSTOM_TARGET(quicktest COMMAND ctest -E "TestR")

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

SET(RTEST scripts/r.sh)
SET(RTESTOUTPUT rtest/regression)

ADD_CUSTOM_TARGET(testall 
  COMMAND ${R_EXECUTABLE} CMD build --no-docs ../../ 
  COMMAND R CMD INSTALL --libs-only qtl_1.13-1.tar.gz 
  COMMAND ctest)

ADD_TEST(TestStandaloneMQM ./scripts/regression_tests.sh)

ADD_TEST(TestR_io ${RTEST} ../../tests test_io.R)
ADD_TEST(TestR_qtl ${RTEST} ../../tests test_qtl.R)
ADD_TEST(TestR_scanone_mr ${RTEST} rtest test_scanone_mr.R)
ADD_TEST(TestR_scanone_mr_cmp ${CMAKE_COMMAND} -E compare_files ${RTESTOUTPUT}/scanone_mr.rnew ${RTESTOUTPUT}/scanone_mr.rtest)
ADD_TEST(TestR_mqm_listeria1 ${RTEST} rtest test_mqm_listeria1.R)
ADD_TEST(TestR_mqm_listeria1_cmp ${CMAKE_COMMAND} -E compare_files ${RTESTOUTPUT}/mqm_listeria1.rnew ${RTESTOUTPUT}/mqm_listeria1.rtest)



