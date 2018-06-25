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

\*---------------------------------------------------------------------------*/

#include "foamVtkFvMeshAdaptor.H"

// VTK includes
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

// Templates (only needed here)
#include "foamVtkFvMeshAdaptorFieldTemplates.C"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
static wordHashSet supportedTypes()
{
    // typedef DimensionedField<Type, volMesh> DimFieldType;
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    wordHashSet types;

    // TODO: types.insert(DimFieldType::typeName);
    types.insert(VolFieldType::typeName);

    return types;
}

} // End namespace Foam


Foam::wordHashSet Foam::vtk::fvMeshAdaptor::knownFields
(
    const wordRes& selectFields
) const
{
    wordHashSet allFields;

    // Quick exit if no volume fields can be converted.
    // This could be refined
    HashTable<wordHashSet> objects = mesh_.classes(selectFields);

    if (objects.empty())
    {
        return allFields;
    }

    wordHashSet types;

    types += supportedTypes<scalar>();
    types += supportedTypes<vector>();
    types += supportedTypes<sphericalTensor>();
    types += supportedTypes<symmTensor>();
    types += supportedTypes<tensor>();

    objects.retain(types);

    forAllConstIters(objects, iter)
    {
        allFields += iter.object();
    }

    return allFields;
}


void Foam::vtk::fvMeshAdaptor::convertVolFields
(
    const wordRes& selectFields
)
{
    // Quick exit if no volume fields can be converted.
    // This could be refined
    HashTable<wordHashSet> objects = mesh_.classes(selectFields);

    if (objects.empty())
    {
        return;
    }

    PtrList<patchInterpolator> interpLst;

    if (interpFields_ && patchIds_.size())
    {
        // NOTE: this would be broken with processor patches,
        // but we don't allow them for the catalyst adaptor anyhow

        // patchIds_ are sorted, so the last one is also the max

        interpLst.setSize(patchIds_.last() + 1);

        for (const label patchId : patchIds_)
        {
            interpLst.set
            (
                patchId,
                new PrimitivePatchInterpolation<primitivePatch>
                (
                    mesh_.boundaryMesh()[patchId]
                )
            );
        }
    }

    convertVolFields<scalar>(interpLst, selectFields);
    convertVolFields<vector>(interpLst, selectFields);
    convertVolFields<sphericalTensor>(interpLst, selectFields);
    convertVolFields<symmTensor>(interpLst, selectFields);
    convertVolFields<tensor>(interpLst, selectFields);

    // TODO
    // convertDimFields<scalar>(interpLst, selectFields);
    // convertDimFields<vector>(interpLst, selectFields);
    // convertDimFields<sphericalTensor>(interpLst, selectFields);
    // convertDimFields<symmTensor>(interpLst, selectFields);
    // convertDimFields<tensor>(interpLst, selectFields);
}


// ************************************************************************* //
