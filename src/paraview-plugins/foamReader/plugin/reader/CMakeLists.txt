#-----------------------------------------------------------------------------
set(classes vtkPVFoamReader)

vtk_module_add_module(FoamReader
  FORCE_STATIC
  CLASSES ${classes})

target_include_directories(FoamReader
    PUBLIC
    ${LIB_SRC}/OpenFOAM/include
    ${LIB_SRC}/OpenFOAM/lnInclude
    ${LIB_SRC}/OSspecific/${WM_OSTYPE}/lnInclude
    ${LIB_SRC}/fileFormats/lnInclude
    ${LIB_SRC}/conversion/lnInclude
    ${LIB_SRC}/finiteArea/lnInclude
    ${LIB_SRC}/finiteVolume/lnInclude
    ${PROJECT_SOURCE_DIR}/../common
    ${PROJECT_SOURCE_DIR}/library
)

set(OPENFOAM_LIBNAMES
    OpenFOAM
    fileFormats
    conversion
    finiteArea
    finiteVolume
)


# With cmake 3.12 and newer, can simply use target_link_directories

# target_link_directories(FoamReader
#     PUBLIC
#     $ENV{FOAM_LIBBIN}
# )

# target_link_libraries(FoamReader
#     PUBLIC
#     vtkPVFoamCommon${foamPvLibQualifier}
#     vtkPVFoamReader${foamPvLibQualifier}
#     ${OPENFOAM_LIBNAMES}
# )


# Older cmake. Do manually

link_directories(
    $ENV{FOAM_LIBBIN}
)

unset(OPENFOAM_LIBS)

foreach(libname ${OPENFOAM_LIBNAMES})
    unset(libpath CACHE)
    find_library(libpath ${libname} PATHS $ENV{FOAM_LIBBIN})
    list(APPEND OPENFOAM_LIBS ${libpath})
endforeach(libname)

target_link_libraries(FoamReader
    PUBLIC
    vtkPVFoamCommon${foamPvLibQualifier}
    vtkPVFoamReader${foamPvLibQualifier}
    ${OPENFOAM_LIBS}
)

#-----------------------------------------------------------------------------
