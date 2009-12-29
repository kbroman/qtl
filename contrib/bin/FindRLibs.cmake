# - Find R libraries - much simplified and works with NIX
#
# This module finds if R is installed and determines where the include files
# and libraries are. It also determines what the name of the library is. This
# code sets the following variables:
#
#   R_INCLUDE_PATH   = path to where R.h is found
#   R_EXECUTABLE     = full path to the R binary
#   R_LIBRARY        = path to shared R library
#   R_BLAS_LIBRARY
#   BIOLIB_R_LIBRARY
#
# They can be hard coded by uncommenting, for example:
#
# SET(R_INCLUDE_PATH "c:/Progra~1/R/R-2.9.1/include")
# SET(R_EXECUTABLE "c:/Progra~1/R/R-2.9.1/bin/R.EXE")
#
# R can be queried with:
#
#   R CMD config --cppflags
#   R CMD config --ldflags
#
# Also pkg-config can be used with --cflags libR etc.

MESSAGE(STATUS,"FindRLibs.cmake")

ASSERT_FOUNDMAP()

FIND_PROGRAM(R_EXECUTABLE R)
IF(WIN32 AND NOT R_EXECUTABLE)
  FIND_PROGRAM(R_EXECUTABLE R
    PATHS "c:/Progra~1/R/R-2.9.1/"
  )
ENDIF()

IF(R_EXECUTABLE)
  GET_FILENAME_COMPONENT(R_BINPATH ${R_EXECUTABLE} PATH)  
  GET_FILENAME_COMPONENT(R_PATH ${R_BINPATH} PATH)  
	SET(R_LIBRARY ${R_PATH}/lib/R/lib/libR.so)
	SET(R_LIKELY_INCLUDE_PATH ${R_PATH}/lib/R/include)
  IF (WIN32)
    SET(R_LIBRARY ${R_BINPATH}/R.DLL)
    SET(R_LIKELY_INCLUDE_PATH ${R_PATH}/include)
  ELSE()
    SET(R_LIBRARY ${R_PATH}/lib/R/lib/libR.so)
    # SET(R_BLAS_LIBRARY ${R_PATH}/lib/R/lib/libR_BLAS.so)
    SET(R_LIKELY_INCLUDE_PATH ${R_PATH}/lib/R/include)
  ENDIF()
ENDIF(R_EXECUTABLE)

# ---- Find R.h and the Rlib.so shared library
SET(R_POSSIBLE_INCLUDE_PATHS
  /usr/share/R/include
)
FIND_PATH(R_INCLUDE_PATH R.h
  ${R_LIKELY_INCLUDE_PATH}
  ${R_POSSIBLE_INCLUDE_PATHS}
)

INCLUDE_DIRECTORIES(${R_INCLUDE_PATH})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(RLibs DEFAULT_MSG R_LIBRARY R_INCLUDE_PATH)

GET_FILENAME_COMPONENT(R_LIBRARY_PATH ${R_LIBRARY} PATH)  

# Locate R_BLAS
FIND_LIBRARY(R_BLAS_LIBRARY
  NAMES Rlbas.dll.a R_BLAS.dll R_BLAS.dylib libR_BLAS.so
  PATHS ${R_LIBRARY_PATH} ${R_BINPATH}
  )

MESSAGE("R_EXECUTABLE=${R_EXECUTABLE}")
MESSAGE("R_INCLUDE_PATH=${R_INCLUDE_PATH}")
MESSAGE("R_LIBRARY=${R_LIBRARY}")
MESSAGE("R_BLAS_LIBRARY=${R_BLAS_LIBRARY}")

if(NOT EXISTS ${R_LIBRARY})
	message(FATAL_ERROR "${R_LIBRARY} was not found (has R been built as shared library libR.so?)")
endif(NOT EXISTS ${R_LIBRARY})

# ---- find the biolib_R Clib module for mapping...
SET(_RLIBNAME ${MAP_projectname}_R)
SET(_RLIBPATH ${MAP_CLIBS_PATH}/${_RLIBNAME})
INCLUDE_DIRECTORIES(${_RLIBPATH}/include)

if(NOT BIOLIB_R_LIBRARY)
  SET(_LIBNAME ${_RLIBNAME}-${MAP_VERSION})
  if(NOT BUILD_LIBS)
    # Look for BIOLIB_R_LIBRARY
    # not building from root, need to find the shared libs
    SET(_LINKLIB lib${_LIBNAME}.so)
    IF(CYGWIN)
      SET(_LINKLIB lib${_LIBNAME}.dll.a)
    ENDIF(CYGWIN)
    IF(APPLE)
      SET(_LINKLIB lib${_LIBNAME}.dylib)
    ENDIF(APPLE)
    message("Looking for ${_LINKLIB} in ${_RLIBPATH}")
    FIND_LIBRARY(R_LIBRARY NAMES ${_LINKLIB} HINTS ${_RLIBPATH}/build ${_RLIBPATH}/src)
    IF(NOT BIOLIB_R_LIBRARY)
      FIND_LIBRARY(BIOLIB_R_LIBRARY NAMES ${_LINKLIB} PATHS ${_RLIBPATH}/build ${_RLIBPATH}/src)
    ENDIF()
    message("Found ${BIOLIB_R_LIBRARY}")
  else()
    SET(BIOLIB_R_LIBRARY ${_LIBNAME})
  endif()
  message("BIOLIB_R_LIBRARY=${BIOLIB_R_LIBRARY}")
endif()
# UNSET(_LIBNAME)
SET(_LIBNAME 'unknown')

MARK_AS_ADVANCED(
  R_INCLUDE_PATH
  R_EXECUTABLE
  R_LIBRARY
  # R_BLAS_LIBRARY
  BIOLIB_R_LIBRARY
  )
