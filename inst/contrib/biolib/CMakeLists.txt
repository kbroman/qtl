# Make file for the R/qtl BioLib module 
#
#   Usage: rm CMakeCache.txt ; cmake . ; make
#
# this cmake file requires the full biolib source tree

cmake_minimum_required(VERSION 2.6)

PROJECT(Biolib_rqtl)

set (M_NAME rqtl)

IF(NOT BUILD_LIBS)
  # Building in a local module: Set the module search path manually
  SET (MAP_ROOT ../../../../..)
  SET (CMAKE_MODULE_PATH ${MAP_ROOT}/tools/cmake-support/modules)
  IF(NOT IS_DIRECTORY ${CMAKE_MODULE_PATH})
    MESSAGE(FATAL_ERROR "Can not find BioLib context")
  ENDIF(NOT IS_DIRECTORY ${CMAKE_MODULE_PATH})
ENDIF(NOT BUILD_LIBS)

FIND_PACKAGE(Map REQUIRED)
FIND_PACKAGE(RLibs REQUIRED)  # require some Rlib functionality at runtime

INCLUDE_DIRECTORIES(../../src)
# INCLUDE_DIRECTORIES(${R_INCLUDE_PATH})

BUILD_CLIB()

ADD_LIBRARY(${LIBNAME} SHARED
../../src/countXO.c
../../src/countXO.h
../../src/discan.c
../../src/discan_covar.c
../../src/discan_covar.h
../../src/discan.h
../../src/effectscan.c
../../src/effectscan.h
../../src/findDupMarkers_notexact.c
../../src/findDupMarkers_notexact.h
../../src/fitqtl_hk.c
../../src/fitqtl_hk.h
../../src/fitqtl_imp.c
../../src/fitqtl_imp.h
../../src/forwsel.c
../../src/forwsel.h
../../src/hmm_4way.c
../../src/hmm_4way.h
../../src/hmm_bc.c
../../src/hmm_bc.h
../../src/hmm_bci.c
../../src/hmm_bci.h
../../src/hmm_f2.c
../../src/hmm_f2.h
../../src/hmm_f2i.c
../../src/hmm_f2i.h
../../src/hmm_main.c
../../src/hmm_main.h
../../src/hmm_ri4self.c
../../src/hmm_ri4self.h
../../src/hmm_ri4sib.c
../../src/hmm_ri4sib.h
../../src/hmm_ri8self.c
../../src/hmm_ri8self.h
../../src/hmm_ri8sib.c
../../src/hmm_ri8sib.h
../../src/info.c
../../src/info.h
../../src/lapackutil.c
../../src/lapackutil.h
../../src/mqmaugment.cpp
../../src/mqmaugment.h
../../src/mqmdatatypes.cpp
../../src/mqmdatatypes.h
../../src/mqmeliminate.cpp
../../src/mqmeliminate.h
../../src/mqmextra.cpp
../../src/mqmextra.h
../../src/mqm.h
../../src/mqmmapqtl.cpp
../../src/mqmmapqtl.h
../../src/mqmmixture.cpp
../../src/mqmmixture.h
../../src/mqmprob.cpp
../../src/mqmprob.h
../../src/mqmregression.cpp
../../src/mqmregression.h
../../src/mqmscan.cpp
../../src/mqmscan.h
../../src/ril48_reorg.c
../../src/ril48_reorg.h
../../src/ripple.c
../../src/ripple.h
../../src/scanone_ehk.c
../../src/scanone_ehk.h
../../src/scanone_em.c
../../src/scanone_em_covar.c
../../src/scanone_em_covar.h
../../src/scanone_em.h
../../src/scanone_hk.c
../../src/scanone_hk.h
../../src/scanone_imp.c
../../src/scanone_imp.h
../../src/scanone_mr.c
../../src/scanone_mr.h
../../src/scanone_np.c
../../src/scanone_np.h
../../src/scantwo_binary_em.c
../../src/scantwo_binary_em.h
../../src/scantwo_em.c
../../src/scantwo_em.h
../../src/scantwo_hk.c
../../src/scantwo_hk.h
../../src/scantwo_imp.c
../../src/scantwo_imp.h
../../src/scantwo_mr.c
../../src/scantwo_mr.h
../../src/simulate.c
../../src/simulate.h
../../src/simulate_ril.c
../../src/simulate_ril.h
../../src/stahl_mf.c
../../src/stahl_mf.h
../../src/standalone.h
../../src/summary_scantwo.c
../../src/summary_scantwo.h
../../src/util.c
../../src/util.h
../../src/vbscan.c
../../src/vbscan.h
)

# ---- The following may be required for OSX
# TARGET_LINK_LIBRARIES(${LIBNAME} ${ZLIB_NAME})


INSTALL_CLIB()
