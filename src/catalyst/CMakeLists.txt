#------------------------------------------------------------------------------
project(catalyst)

cmake_minimum_required(VERSION 3.8)
cmake_policy(SET CMP0002 NEW)   # For cmake >= 3
cmake_policy(SET CMP0057 NEW)   # For cmake >= 3.10.2
set(CMAKE_BUILD_TYPE Release)

# Fail if not building out-of-source
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds disallowed. Use a separate build directory")
endif()

#------------------------------------------------------------------------------
# Simple discovery and sanity checks

if (EXISTS "$ENV{ParaView_DIR}")
    find_package(ParaView REQUIRED)
else()
    message(FATAL_ERROR "ParaView_DIR not set")
endif()

# Define pv_api as MAJ.MIN
if (${ParaView_VERSION} VERSION_GREATER_EQUAL 5.7)
    string(REGEX MATCH [0-9]+[.][0-9]+ pv_api ${ParaView_VERSION})
    set(pv_version "${ParaView_VERSION}")
else()
    set(pv_api "${PARAVIEW_VERSION_MAJOR}.${ParaView_VERSION_MINOR}")
    set(pv_version "${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}.${PARAVIEW_VERSION_PATCH}")
endif()

if (${ParaView_VERSION} VERSION_GREATER_EQUAL 5.7)
    find_package(ParaView REQUIRED
        COMPONENTS PythonCatalyst
        OPTIONAL_COMPONENTS VTK::ParallelMPI)
else()
    find_package(ParaView REQUIRED
        COMPONENTS vtkPVPythonCatalyst
        OPTIONAL_COMPONENTS vtkParallelMPI)
    ## include(${PARAVIEW_USE_FILE})
endif()

include(GNUInstallDirs)

#-----------------------------------------------------------------------------

message("================")
message("Using OpenFOAM = $ENV{WM_PROJECT_DIR}")
message("Using ParaView = $ENV{ParaView_DIR}")
if (NOT PARAVIEW_USE_MPI)
    message(WARNING "==== Recommended to build using ParaView Catalyst with MPI ====")
endif()

include(CMakeLists-OpenFOAM.txt)

# Default is PROJECT target, with supporting libraries from FOAM_LIBBIN
# Possible temporary override with 'STAGED_INSTALL_PREFIX' later
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    get_filename_component(FOAM_PLATFORMS ${FOAM_LIBBIN} DIRECTORY)
    set(CMAKE_INSTALL_PREFIX "${FOAM_PLATFORMS}"
        CACHE PATH "Install path prefix" FORCE)
endif()
message("Install prefix = ${CMAKE_INSTALL_PREFIX}")


message("================")
include(CMakeLists-Project.txt)

#-----------------------------------------------------------------------------
