#------------------------------------------------------------------------------
project(runTimePostProcessing)

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

set(ParaView_VERSION 0)
unset(VTK_VERSION)
unset(VTK_VERSION_MAJOR)

if (EXISTS "$ENV{VTK_DIR}")
    message("Building with VTK from $ENV{VTK_DIR}")
    find_package(VTK REQUIRED)

elseif (EXISTS "$ENV{ParaView_DIR}")
    message("Building with Paraview from $ENV{ParaView_DIR}")

    find_package(ParaView REQUIRED)
    if (${ParaView_VERSION} VERSION_GREATER_EQUAL 5.7)
        find_package(ParaView REQUIRED
            OPTIONAL_COMPONENTS VTK::ParallelMPI)
    endif()
else()
    message(FATAL_ERROR "VTK not found using VTK_DIR or ParaView_DIR")
endif()

if (VTK_USE_FILE)
    message("VTK from ${VTK_USE_FILE}")
    include(${VTK_USE_FILE})
endif()

include(GNUInstallDirs)

#-----------------------------------------------------------------------------
# Handle name changes (eg, VTK_MAJOR_VERSION to VTK_VERSION_MAJOR etc.)

if (VTK_MAJOR_VERSION AND NOT VTK_VERSION_MAJOR)
    message("Found older VTK version naming")
    set(VTK_VERSION_MAJOR ${VTK_MAJOR_VERSION})
    set(VTK_VERSION_MINOR ${VTK_MINOR_VERSION})
    set(VTK_VERSION_PATCH ${VTK_BUILD_VERSION})

    if (NOT VTK_VERSION)
        set(
            VTK_VERSION
            "${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}.${VTK_VERSION_PATCH}"
        )
        message("Synthesized VTK version: " ${VTK_VERSION})
    endif()
endif()

if (${ParaView_VERSION} GREATER 0)
    # Newer VTK - try to load optional MPI?
    if (${VTK_VERSION} GREATER 8.2)
        find_package(VTK REQUIRED OPTIONAL_COMPONENTS ParallelMPI)
    endif()
endif()

#-----------------------------------------------------------------------------

message("================")
message("Using OpenFOAM = $ENV{WM_PROJECT_DIR}")
if (${ParaView_VERSION} GREATER 0)
message("Using ParaView = ${ParaView_VERSION} (VTK ${VTK_VERSION})")
else()
message("Using VTK ${VTK_VERSION}")
endif()

include(CMakeLists-OpenFOAM.txt)
include(CMakeLists-library.txt)

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
