#------------------------------------------------------------------------------
project(ParaFoamBlockReader)

cmake_minimum_required(VERSION 3.8)
cmake_policy(SET CMP0002 NEW)   # For cmake >= 3
cmake_policy(SET CMP0057 NEW)   # For cmake >= 3.10.2
set(CMAKE_BUILD_TYPE Release)

# Fail if not building out-of-source
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds disallowed. Use a separate build directory")
endif()

#-----------------------------------------------------------------------------
# Simple discovery and sanity checks

find_package(ParaView REQUIRED)

# Qualify OpenFOAM/ParaView libraries with -pvMAJ.MIN
# The pv_api can also be used for paraview-MAJ.MIN
if (DEFINED ENV{PARAVIEW_API})
    set(pv_api "$ENV{PARAVIEW_API}")
else()
    if (${ParaView_VERSION} VERSION_GREATER_EQUAL 5.7)
        string(REGEX MATCH [0-9]+[.][0-9]+ pv_api ${ParaView_VERSION})
        message("The PARAVIEW_API env was unset, "
            "extracted value from version="
            "${ParaView_VERSION} is ${pv_api}")
    else()
        set(pv_api "${PARAVIEW_VERSION_MAJOR}.${ParaView_VERSION_MINOR}")
        message("The PARAVIEW_API env was unset, "
            "assembled from maj/min="
            "${ParaView_VERSION_MAJOR}/${ParaView_VERSION_MINOR} is ${pv_api}")
    endif()
endif()

# Pass through the OpenFOAM API value to the reader
set(foam_api "$ENV{FOAM_API}")
if (NOT "${foam_api}")
    set(foam_api "2.0")
endif()

set(foamPvLibQualifier "-pv${pv_api}")

include(GNUInstallDirs)

#-----------------------------------------------------------------------------

message("================")
message("Using ParaView = $ENV{ParaView_DIR}")
message("OpenFOAM api   = ${foam_api}, libs = ${foamPvLibQualifier}")

#DEBUG get_cmake_property(_variableNames VARIABLES)
#DEBUG foreach (_variableName ${_variableNames})
#DEBUG     message(STATUS "${_variableName}=${${_variableName}}")
#DEBUG endforeach()

#-----------------------------------------------------------------------------

# Default is PROJECT target, with supporting libraries from FOAM_LIBBIN
# Possible temporary override with 'STAGED_INSTALL_PREFIX' later
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    get_filename_component(FOAM_PLATFORMS ${FOAM_LIBBIN} DIRECTORY)
    set(CMAKE_INSTALL_PREFIX "${FOAM_PLATFORMS}"
        CACHE PATH "Install path prefix"
        FORCE)
    set(SUPPORTING_LIBBIN "${FOAM_LIBBIN}")
else()
    set(SUPPORTING_LIBBIN "${CMAKE_INSTALL_PREFIX}/lib")
endif()
message("Install prefix = ${CMAKE_INSTALL_PREFIX}")

include(CMakeLists-OpenFOAM.txt)
include(CMakeLists-library.txt)

if (${ParaView_VERSION} VERSION_GREATER_EQUAL 5.7)
  message("================")
  include(CMakeLists-Project57.txt)
else()
  message("================")
  include(${PARAVIEW_USE_FILE})
  include(CMakeLists-Project.txt)
endif()

#-----------------------------------------------------------------------------
