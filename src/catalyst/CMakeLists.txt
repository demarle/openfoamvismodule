cmake_minimum_required(VERSION 2.8)

project(catalyst)

# Set policy for CMP0002 needed for cmake > 3
cmake_policy(SET CMP0002 OLD)

if (EXISTS "$ENV{ParaView_DIR}")
    message("Building with Paraview from $ENV{ParaView_DIR}")
    find_package(ParaView REQUIRED COMPONENTS vtkPVPythonCatalyst)
    if (NOT PARAVIEW_USE_MPI)
        message(WARNING "==== Recommended to build using ParaView Catalyst with MPI ====")
    endif()
else()
    message (FATAL_ERROR "ParaView_DIR not set")
endif()

include(CMakeLists-Common.txt)


#-----------------------------------------------------------------------------
