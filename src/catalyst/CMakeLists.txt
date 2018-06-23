#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)

# Fail immediately if not using an out-of-source build
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds are not supported.  Please create a build directory "
        "separate from the source directory")
endif()

project(catalyst)

# Set policy for CMP0002 needed for cmake > 3
cmake_policy(SET CMP0002 OLD)


#------------------------------------------------------------------------------
# Simple sanity checks

if (NOT EXISTS "$ENV{WM_PROJECT_DIR}")
    message(FATAL_ERROR "WM_PROJECT_DIR not set")
endif()

if (EXISTS "$ENV{ParaView_DIR}")
    find_package(ParaView REQUIRED COMPONENTS vtkPVPythonCatalyst)
else()
    message(FATAL_ERROR "ParaView_DIR not set")
endif()


#------------------------------------------------------------------------------
# Installation locations (not much)
# - default install into FOAM_LIBBIN

if ((NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY) OR
    ("${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" STREQUAL "default"))
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY $ENV{FOAM_LIBBIN})
endif()


#-----------------------------------------------------------------------------
# All messages

message("================")
message("Using OpenFOAM = $ENV{WM_PROJECT_DIR}")
message("Using ParaView = $ENV{ParaView_DIR}")
if (NOT PARAVIEW_USE_MPI)
    message(WARNING "==== Recommended to build using ParaView Catalyst with MPI ====")
endif()
message("install libdir = ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message("================")


#-----------------------------------------------------------------------------

include(CMakeLists-Common.txt)

#-----------------------------------------------------------------------------
