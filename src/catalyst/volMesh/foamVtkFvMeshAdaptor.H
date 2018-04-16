/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::vtk::fvMeshAdaptor

Description
    The backend for the catalystFvMesh function object.

    The output is a multi-block dataset with two blocks corresponding
    to the internal (volume) mesh (block 0) and the patches as (block 1).
    These two blocks correspond simultaneously to the catalyst channels
    "mesh" and "patches", respectively.

    The patches are further divided into sub-blocks. The lowest data blocks
    are multi-piece datasets (UnstructuredGrid and PolyData for internal mesh
    and patches, respectively) with each piece corresponding to its MPI rank.

SourceFiles
    foamVtkFvMeshAdaptor.C
    foamVtkFvMeshAdaptorGeom.C
    foamVtkFvMeshAdaptorGeomVtu.C
    foamVtkFvMeshAdaptorFields.C
    foamVtkFvMeshAdaptorFieldTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkFvMeshAdaptor_H
#define foamVtkFvMeshAdaptor_H

#include "className.H"
#include "fileName.H"
#include "stringList.H"
#include "wordList.H"
#include "Enum.H"
#include "primitivePatch.H"
#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"
#include "foamVtkTools.H"
#include "foamVtkMeshMaps.H"

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMultiBlockDataSet.h>

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

class vtkCellArray;
class vtkDataSet;
class vtkFloatArray;
class vtkIndent;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;
template<class Type> class Field;

namespace vtk
{

/*---------------------------------------------------------------------------*\
                     Class vtk::fvMeshAdaptor Declaration
\*---------------------------------------------------------------------------*/

class fvMeshAdaptor
:
    protected vtk::Tools
{
public:

        //- Public Data Members

        //- The Catalyst output channels
        enum channel
        {
            MESH = 0x1,
            PATCHES = 0x2,
            INPUT = 0x3
        };

        static const Enum<channel> channelNames;

        //- Name for internal mesh ("internal")
        static const word internalName;


private:

    // Convenience typedefs
    typedef PrimitivePatchInterpolation<primitivePatch> patchInterpolator;

        //- Bookkeeping for vtkPolyData
        struct foamVtpData
        :
            public vtk::Caching<vtkPolyData>,
            public foamVtkMeshMaps
        {};

        //- Bookkeeping for vtkUnstructuredGrid
        struct foamVtuData
        :
            public vtk::Caching<vtkUnstructuredGrid>,
            public foamVtkMeshMaps
        {

            //- The vtk points for the mesh (and decomposition)
            vtkSmartPointer<vtkPoints> points
            (
                const fvMesh& mesh
            ) const;

            //- The vtk points for the mesh (and decomposition)
            //- using the provided pointMap
            vtkSmartPointer<vtkPoints> points
            (
                const fvMesh& mesh,
                const labelUList& pointMap
            ) const;

            //- Internal mesh as vtkUnstructuredGrid
            vtkSmartPointer<vtkUnstructuredGrid> internal
            (
                const fvMesh& mesh,
                const bool decompPoly
            );
        };


    // Private Data

        //- OpenFOAM mesh
        const fvMesh& mesh_;

        //- Selected output channels
        unsigned channels_;

        //- Previous/current decomposition request
        bool decomposePoly_;

        //- Track changes in mesh geometry
        enum polyMesh::readUpdateState meshState_;

        //- Any information for 2D (VTP) geometries
        HashTable<foamVtpData, string> cachedVtp_;

        //- Cell maps and other information for 3D (VTU) geometries
        HashTable<foamVtuData, string> cachedVtu_;


    // Mesh Conversion

        //- Convert internal
        void convertGeometryInternal();

        //- Convert patches
        //  Patches may use additionalIds (cached data) with the patch Ids.
        //  There will be several for groups, but only one for regular patches.
        void convertGeometryPatches();

        //- Ghost cells to affect the visibility of the geometry.
        //  Currently only used for overset.
        void applyGhosting();


    // Field Conversion

        //- Convert specified volume fields
        void convertVolFields(const wordRes& selectFields);

        //- Volume field
        template<class Type>
        vtkSmartPointer<vtkFloatArray> convertVolFieldToVTK
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const foamVtuData& vtuData
        ) const;


        //- Volume field - all types
        template<class Type>
        void convertVolField
        (
            const PtrList<patchInterpolator>& patchInterpList,
            const GeometricField<Type, fvPatchField, volMesh>& fld
        );

        //- Volume fields - all types
        template<class Type>
        void convertVolFields
        (
            const PtrList<patchInterpolator>& patchInterpList,
            const wordRes& selectFields
        );

        //- Volume internal fields (DimensionedField)- all types
        template<class Type>
        void convertDimFields
        (
            const PtrList<patchInterpolator>& patchInterpList,
            const wordRes& selectFields
        );

        //- Volume field - all selected parts
        template<class Type>
        void convertVolFieldInternal
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            autoPtr<GeometricField<Type, pointPatchField, pointMesh>>& ptfPtr
        );

        //- Point field
        template<class Type>
        vtkSmartPointer<vtkFloatArray> convertPointField
        (
            const GeometricField<Type, pointPatchField, pointMesh>& pfld,
            const GeometricField<Type, fvPatchField, volMesh>& vfld,
            const foamVtuData& vtuData
        );


        //- Update geometry and fields
        void updateContent(const wordRes& selectFields);


    // Constructors

        //- No copy construct
        fvMeshAdaptor(const fvMeshAdaptor&) = delete;

        //- No copy assignment
        void operator=(const fvMeshAdaptor&) = delete;


public:

    //- Static Data Members

        ClassName("vtk::fvMeshAdaptor");


    // Constructors

        //- Construct from components
        fvMeshAdaptor(const fvMesh& mesh);


    //- Destructor
    ~fvMeshAdaptor() = default;


    // Member Functions

        //- Define the output channels by name
        void channels(const wordList& chanNames);

        //- Define the output channels by value
        void channels(unsigned chanIds);

        //- Return the selected output channel ids
        label channels() const;

        //- True if MESH channel is being used
        bool usingVolume() const;

        //- True if PATCHES channel is being used
        bool usingPatches() const;

        //- Number of non-processor patches, when the PATCHES channel is being
        //- used - otherwise 0.
        label nPatches() const;

        //- Return the names of known (supported) fields
        wordHashSet knownFields(const wordRes& selectFields) const;

        void updateState(polyMesh::readUpdateState state);

        //- The output is a pair (MESH/PATCHES channels) of vtk meshes with
        //- point/cell data.
        vtkSmartPointer<vtkMultiBlockDataSet> output
        (
            const wordRes& selectFields
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
