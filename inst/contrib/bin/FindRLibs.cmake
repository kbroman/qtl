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
IF(WIN32) 
	SET(R_INCLUDE_PATH "C:/Progra~1/R/R-3.1.0/include")
	SET(R_EXECUTABLE "C:/Progra~1/R/R-3.1.0/bin/i386/R.exe")
	SET(R_LIBRARY "C:/Progra~1/R/R-3.1.0/bin/i386/R.dll")
	SET(R_BLAS_LIBRARY "C:/Progra~1/R/R-3.1.0/bin/i386/Rblas.dll")
ENDIF()
#
# R can be queried with:
#
#   R CMD config --cppflags - returns '-I/usr/share/R/include -I/usr/share/R/include'
#   R CMD config --ldflags  - returns '-L/usr/lib/R/lib -lR'
#
# Also pkg-config can be used with --cflags libR etc.

MESSAGE("FindRLibs.cmake")

ASSERT_FOUNDMAP()

IF(APPLE)
  message("OSX environment")
  FIND_PROGRAM(R_EXECUTABLE R
    NO_DEFAULT_PATH
    PATHS 
      /Library/Frameworks/R.framework/Resources/bin
      /usr/bin
      /usr/local/bin
  )
ENDIF()

# If the R path is in the cache the following will skip...
FIND_PROGRAM(R_EXECUTABLE R
  PATHS      # maybe use HINTS
    /usr/bin
    /usr/local/bin
)

IF(WIN32 AND NOT R_EXECUTABLE)
  FIND_PROGRAM(R_EXECUTABLE R
    PATHS "c:/Progra~1/R/R-3.1.0"
  )
ENDIF()

IF(R_EXECUTABLE)
  message("R_EXECUTABLE=${R_EXECUTABLE}")
  GET_FILENAME_COMPONENT(R_BINPATH ${R_EXECUTABLE} PATH)  
  GET_FILENAME_COMPONENT(R_PATH ${R_BINPATH} PATH)
  # Get information from R itself
  # Fetch the library paths
  IF(NOT WIN32)
    EXECUTE_PROCESS(COMMAND ${R_EXECUTABLE} CMD config --ldflags OUTPUT_VARIABLE _LIBS)
  ELSE()
    SET(_LIBS "-LC:/PROGRA~1/R/R-3.1.0/bin -lR")
  ENDIF()
  message("LIBS=${_LIBS}")
  if (APPLE)
    # returns -F/Library/Frameworks/R.framework/.. -framework R
    STRING(REGEX REPLACE "-F([^ ]+) -framework R" "\\1" R_EXE_LIB_PATHS "${_LIBS}")
  else()
    STRING(REGEX REPLACE "-L([^ ]+)" "\\1" R_EXE_LIB_PATHS "${_LIBS}")
  endif()
  separate_arguments(R_EXE_LIB_PATHS)
  message("R_EXE_LIB_PATHS=${R_EXE_LIB_PATHS}")

  # Fetch the include paths
  IF(NOT WIN32)  
    EXECUTE_PROCESS(COMMAND ${R_EXECUTABLE} CMD config --cppflags OUTPUT_VARIABLE _INCLUDES)
  ELSE()
    SET(_INCLUDES "-IC:/PROGRA~1/R/R-3.1.0/include")
  ENDIF()
  message("INCLUDES=${_INCLUDES}")
  STRING(REGEX REPLACE "-I([^ ]+)" "\\1" R_EXE_INCLUDE_PATHS "${_INCLUDES}")
  separate_arguments(R_EXE_INCLUDE_PATHS)
  message("R_EXE_INCLUDE_PATHS=${R_EXE_INCLUDE_PATHS}")

  # Locate libraries e.g. Linux /usr/lib/R/lib/libR.so /usr/bin/R
  FIND_LIBRARY(R_LIBRARY
    NAMES libR.so R.DLL libR.dylib
    PATHS 
      ${R_EXE_LIB_PATHS}
      ${R_PATH}
      ${R_BINPATH}
      /Library/Frameworks/R.framework/Resources/lib
  )
	SET(R_LIKELY_INCLUDE_PATH ${R_PATH}/lib/R/include)
  IF (WIN32)
    SET(R_LIKELY_INCLUDE_PATH ${R_PATH}/include)
  ENDIF()
ENDIF(R_EXECUTABLE)

GET_FILENAME_COMPONENT(R_LIBRARY_PATH ${R_LIBRARY} PATH)  

# ---- Find R.h include file(s)
FIND_PATH(R_INCLUDE_PATH R.h
  ${R_EXE_INCLUDE_PATHS}
  ${R_LIKELY_INCLUDE_PATH}
)
INCLUDE_DIRECTORIES(${R_INCLUDE_PATH})

FIND_LIBRARY(R_BLAS_LIBRARY
  NAMES Rlbas.dll.a R_BLAS.dll libRblas.dylib libR_BLAS.so libblas.so
  PATHS 
      ${R_LIBRARY_PATH} 
      ${R_EXE_LIB_PATHS}
      ${R_PATH}
      ${R_BINPATH}
      /Library/Frameworks/R.framework/Resources/lib
  )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(RLibs DEFAULT_MSG R_LIBRARY R_INCLUDE_PATH)

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
    FIND_LIBRARY(R_LIBRARY NAMES ${_LINKLIB} 
      HINTS ${_RLIBPATH}/build ${_RLIBPATH}/src)
    IF(NOT BIOLIB_R_LIBRARY)
      FIND_LIBRARY(BIOLIB_R_LIBRARY NAMES ${_LINKLIB} 
        PATHS ${_RLIBPATH}/build ${_RLIBPATH}/src)
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
  BIOLIB_R_LIBRARY
  )
