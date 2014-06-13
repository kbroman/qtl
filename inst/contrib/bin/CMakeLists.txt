cmake_minimum_required(VERSION 2.6)

IF (NOT CMAKE_MODULE_PATH)
  set(CMAKE_MODULE_PATH ".")
ENDIF()

SET(CMAKE_BUILD_TYPE Debug)
IF(NOT WIN32)  
  EXECUTE_PROCESS(COMMAND grep Version ../../DESCRIPTION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE _VERSION)
  STRING(REGEX REPLACE "Version: " "" RQTL_VERSION "${_VERSION}")
ELSE()
  # Cannot assume we have grep on Windows - should be replaced with
  # CMake logic (FIXME)
  SET(RQTL_VERSION 1.15-3)
ENDIF()
MESSAGE("RQTL_VERSION=${RQTL_VERSION}")

IF ( CMAKE_BUILD_TYPE MATCHES Debug )
  ADD_DEFINITIONS(-Wall -D_DEBUG)
ENDIF ( CMAKE_BUILD_TYPE MATCHES Debug )

# CMake testing
ENABLE_TESTING()
ADD_CUSTOM_TARGET(quicktest COMMAND ctest -E "TestR")


if(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
  set(BITS64 TRUE)
endif()

# ---- The following MACRO is just to keep the FindRLibs.cmake happy
MACRO(ASSERT_FOUNDMAP)
  SET(BUILD_LIBS TRUE)
  SET(BIOLIB_R_LIBRARY UNUSED)
ENDMACRO()

FIND_PACKAGE(RLibs)

include_directories (../../../src)

# Compiler options
add_definitions(-DSTANDALONE)
# ADD_DEFINITIONS(-pedantic)
ADD_DEFINITIONS(-DENABLE_C99_MACROS)


IF(NOT WIN32)
  add_definitions(-pg)
  set(LINK_FLAGS "-pg")
ENDIF()

link_libraries(${R_LIBRARY} ${R_BLAS_LIBRARY})

add_executable (mqm
../../../src/simulate.c
../../../src/simulate.h
../../../src/mqmaugment.cpp    
../../../src/mqmeliminate.h  
../../../src/mqmmapqtl.h     
../../../src/mqmregression.cpp
../../../src/mqmaugment.h      
../../../src/mqmmixture.cpp  
../../../src/mqmregression.h
../../../src/mqmdatatypes.cpp 
../../../src/mqmmixture.h    
../../../src/mqmscan.cpp
../../../src/mqmdatatypes.h    
../../../src/mqm.h           
../../../src/mqmprob.cpp     
../../../src/mqmscan.h
../../../src/mqmeliminate.cpp  
../../../src/mqmmapqtl.cpp   
../../../src/mqmprob.h       
../../../src/standalone.h
../../../src/util.h
../../../src/util.c
mqmmain.cpp
)

add_executable (mqmdebugout
mqmdebugout.cpp
)

SET(RTEST scripts/r.sh)
SET(RTESTOUTPUT rtest/regression)

# Test R version (updating only the C libraries)
ADD_CUSTOM_TARGET(testR_Clibs
  COMMAND ${R_EXECUTABLE} CMD build --no-docs ../../ 
  COMMAND ${R_EXECUTABLE} CMD INSTALL --libs-only qtl_${RQTL_VERSION}.tar.gz 
  COMMAND ${CMAKE_COMMAND} -E remove qtl_${RQTL_VERSION}.tar.gz
  COMMAND ctest -R TestR)

# Test R version only
ADD_CUSTOM_TARGET(testR
  COMMAND ${R_EXECUTABLE} CMD build --no-docs ../../ 
  COMMAND ${R_EXECUTABLE} CMD INSTALL qtl_${RQTL_VERSION}.tar.gz 
  COMMAND ${CMAKE_COMMAND} -E remove qtl_${RQTL_VERSION}.tar.gz
  COMMAND ctest -R TestR)
  
# All tests, including R
ADD_CUSTOM_TARGET(testall
  COMMAND ${R_EXECUTABLE} CMD build --no-docs ../../ 
  COMMAND ${R_EXECUTABLE} CMD INSTALL qtl_${RQTL_VERSION}.tar.gz 
  COMMAND ${CMAKE_COMMAND} -E remove qtl_${RQTL_VERSION}.tar.gz
  COMMAND ctest)

IF(WIN32)
  SET(MQMEXE "mqm.exe")
  SET(MQMDEBUGOUT "mqmdebugout.exe")
ELSE()
  SET(MQMEXE "./mqm")
  SET(MQMDEBUGOUT "./mqmdebugout")
ENDIF()


# Tracetest
ADD_CUSTOM_TARGET(tracetest
  COMMAND ${MQMEXE} -o T12trace.log -d1 -ptest/std/phenotypes1.txt -gtest/std/genotypes1.txt -mtest/std/markers1.txt -stest/std/settings1.txt -ctest/t12/cofactors.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --minprob=1 
)


# debugout output
ADD_TEST(debugout_dnorm ${MQMDEBUGOUT} -d -o test/regression/debugout_dnorm.new)
ADD_TEST(debugout_dnorm_cmp ${CMAKE_COMMAND} -E compare_files test/regression/debugout_dnorm.new test/regression/debugout_dnorm.txt)

ADD_TEST(debugout_pbeta ${MQMDEBUGOUT} -p -o test/regression/debugout_pbeta.new)
ADD_TEST(debugout_pbeta_cmp ${CMAKE_COMMAND} -E compare_files test/regression/debugout_pbeta.new test/regression/debugout_pbeta.txt)

# Older test for comparison with (older test) MQM_test0.txt
# to only run this test use 'ctest -R test0'
ADD_TEST(TestMQMF2_test0 ${MQMEXE} -v -ptest/std/phenotypes1.txt -gtest/std/genotypes1.txt -mtest/std/markers1.txt -stest/std/settings1.txt --smin=-20 --smax=220 --sstep=5 --alpha=0.05 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --minprob=1 -o test/regression/t11out-test0.txt.rnew)
ADD_TEST(TestMQMF2_test0_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t11out-test0.txt.rnew test/regression/t11out-test0.txt)

# Run T11
ADD_TEST(TestMQMF2_T11 ${MQMEXE} -v -ptest/std/phenotypes1.txt -gtest/std/genotypes1.txt -mtest/std/markers1.txt -stest/std/settings1.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --minprob=1 -o test/regression/t11out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMF2_T11_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t11out.txt.rnew test/regression/t11out.txt)

# Run T12 Dataset 1 F2
ADD_TEST(TestMQMF2_T12 ${MQMEXE} -v -ptest/std/phenotypes1.txt -gtest/std/genotypes1.txt -mtest/std/markers1.txt -stest/std/settings1.txt -ctest/t12/cofactors.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --minprob=1 -o test/regression/t12out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMF2_T12_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t12out.txt.rnew test/regression/t12out.txt)

# Run T13 
ADD_TEST(TestMQMF2_T13 ${MQMEXE} -v -ptest/std/phenotypes1.txt -gtest/std/genotypes1.txt -mtest/std/markers1.txt -stest/std/settings1.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=y --maugment=10000 --miaugment=250 --minprob=1 -o test/regression/t13out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMF2_T13_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t13out.txt.rnew test/regression/t13out.txt)

# Run T21 Dataset 2 BC=~ Hyper (genotypes NA -> 0)
ADD_TEST(TestMQMBC_T21 ${MQMEXE} -v -ptest/std/phenotypes2.txt -gtest/std/genotypes2.txt -mtest/std/markers2.txt -stest/std/settings2.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --minprob=1 -o test/regression/t21out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMBC_T21_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t21out.txt.rnew test/regression/t21out.txt)

# Run T22
ADD_TEST(TestMQMBC_T22 ${MQMEXE} -v -ptest/std/phenotypes2.txt -gtest/std/genotypes2.txt -mtest/std/markers2.txt -stest/std/settings2.txt -ctest/t22/cofactors.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --minprob=1 -o test/regression/t22out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMBC_T22_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t22out.txt.rnew test/regression/t22out.txt)

# Run T23
#crash on est.map
ADD_TEST(TestMQMBC_T23 ${MQMEXE} -v -ptest/std/phenotypes2.txt -gtest/std/genotypes2.txt -mtest/std/markers2.txt -stest/std/settings2.txt -ctest/t23/cofactors.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=y --maugment=10000 --miaugment=250 --minprob=1 -o test/regression/t23out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMNC_T23_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t23out.txt.rnew test/regression/t23out.txt)

# Run T24
#crash due to augmentation
ADD_TEST(TestMQMBC_T24 ${MQMEXE} -v -ptest/std/phenotypes2.txt -gtest/std/genotypes2m.txt -mtest/std/markers2.txt -stest/std/settings2.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --minprob=1 -o test/regression/t24out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMBC_T24_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t24out.txt.rnew test/regression/t24out.txt)

# Run T25
ADD_TEST(TestMQMBC_T25 ${MQMEXE} -v -ptest/std/phenotypes2.txt -gtest/std/genotypes2m.txt -mtest/std/markers2.txt -stest/std/settings2.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --minprob=0.5 -o test/regression/t25out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMBC_T25_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t25out.txt.rnew test/regression/t25out.txt)

# Run T31 Dataset 3 F2=~ Listeria
ADD_TEST(TestMQMF2_T31 ${MQMEXE} -v -ptest/std/phenotypes3.txt -gtest/std/genotypes3.txt -mtest/std/markers3.txt -stest/std/settings3.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --minprob=0.5 -o test/regression/t31out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMF2_T31_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t31out.txt.rnew test/regression/t31out.txt)

# Run T32
#if(BITS64)
#  MESSAGE("Skip T32 on 64-bits")
#else()
  ADD_TEST(TestMQMF2_T32 ${MQMEXE} -v -ptest/std/phenotypes3.txt -gtest/std/genotypes3m.txt -mtest/std/markers3.txt -stest/std/settings3.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --minprob=0.5 -o test/regression/t32out.txt.rnew)
  # Immediately compare outputs against older results
  ADD_TEST(TestMQMF2_T32_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t32out.txt.rnew test/regression/t32out.txt)
#endif()

# Run T33
ADD_TEST(TestMQMF2_T33 ${MQMEXE} -v -ptest/std/phenotypes3.txt -gtest/std/genotypes3.txt -mtest/std/markers3.txt -stest/std/settings3.txt -ctest/t33/cofactors.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=10000 --miaugment=250 --minprob=0.5 -o test/regression/t33out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMF2_T33_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t33out.txt.rnew test/regression/t33out.txt)

# Run T32
ADD_TEST(TestMQMF2_T34 ${MQMEXE} -v -ptest/std/phenotypes3.txt -gtest/std/genotypes3m.txt -mtest/std/markers3.txt -stest/std/settings3.txt --smin=0 --smax=200 --sstep=2 --alpha=0.02 --window=10 --maxiter=1000 --estmap=n --maugment=120 --miaugment=1 --minprob=1 -o test/regression/t34out.txt.rnew)
# Immediately compare outputs against older results
ADD_TEST(TestMQMF2_T34_cmp ${CMAKE_COMMAND} -E compare_files test/regression/t34out.txt.rnew test/regression/t34out.txt)

IF(TEST_R)
  ADD_TEST(TestR_io ${RTEST} ../../tests test_io.R)
  ADD_TEST(TestR_qtl ${RTEST} ../../tests test_qtl.R)
  ADD_TEST(TestR_scanone_mr ${RTEST} rtest test_scanone_mr.R)
  ADD_TEST(TestR_scanone_mr_cmp ${CMAKE_COMMAND} -E compare_files ${RTESTOUTPUT}/scanone_mr.rnew ${RTESTOUTPUT}/scanone_mr.rtest)
  ADD_TEST(TestR_mqm_listeria1 ${RTEST} rtest test_mqm_listeria1.R)
  if(NOT BITS64)
    ADD_TEST(TestR_mqm_listeria1_cmp ${CMAKE_COMMAND} -E compare_files ${RTESTOUTPUT}/mqm_listeria1.rnew ${RTESTOUTPUT}/mqm_listeria1.rtest)
  endif()
  ADD_TEST(TestR_augmentation ${RTEST} rtest test_augmentation.R)
ENDIF()


