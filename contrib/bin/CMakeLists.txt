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

IF(NOT WIN32)
  add_definitions(-pg)
  set(LINK_FLAGS "-pg")
  link_libraries(${R_LIBRARY})
ELSE()
  link_libraries(${R_LIBRARY} ${RBLAS_LIBRARY})
ENDIF()


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

IF(WIN32)
  SET(MQMEXE "mqm.exe")
ELSE()
  SET(MQMEXE "./mqm")
ENDIF()

# Older test for comparison with (older test) MQM_test0.txt
ADD_TEST(TestMQMF2_test0 ${MQMEXE} -v -ptest/std/phenotypes1.txt -gtest/std/genotypes1.txt -mtest/std/markers1.txt -stest/std/settings1.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --neglect=1 -o test/regression/t11out-test0.txt.rnew -ctest/t11/cofactors.txt)
ADD_TEST(TestMQMF2_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t11out-test0.txt.rnew test/regression/t11out-test0.txt)

# Run T11
ADD_TEST(TestMQMF2_T11 ${MQMEXE} -v -ptest/std/phenotypes1.txt -gtest/std/genotypes1.txt -mtest/std/markers1.txt -stest/std/settings1.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --neglect=1 -o test/regression/t11out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMF2_T11_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t11out.txt.rnew test/regression/t11out.txt)


# ADD_TEST(DannyTestMQM ./regressiontests.bat) - move into above...

IF(TEST_R)
  ADD_TEST(TestR_io ${RTEST} ../../tests test_io.R)
  ADD_TEST(TestR_qtl ${RTEST} ../../tests test_qtl.R)
  ADD_TEST(TestR_scanone_mr ${RTEST} rtest test_scanone_mr.R)
  ADD_TEST(TestR_scanone_mr_cmp ${CMAKE_COMMAND} -E compare_files ${RTESTOUTPUT}/scanone_mr.rnew ${RTESTOUTPUT}/scanone_mr.rtest)
  ADD_TEST(TestR_mqm_listeria1 ${RTEST} rtest test_mqm_listeria1.R)
  ADD_TEST(TestR_mqm_listeria1_cmp ${CMAKE_COMMAND} -E compare_files ${RTESTOUTPUT}/mqm_listeria1.rnew ${RTESTOUTPUT}/mqm_listeria1.rtest)
ENDIF()


