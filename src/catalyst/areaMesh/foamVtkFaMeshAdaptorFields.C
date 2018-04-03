/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "foamVtkFaMeshAdaptor.H"
#include "faMesh.H"
#include "areaFields.H"

// VTK includes
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

// Templates (only needed here)
#include "foamVtkFaMeshAdaptorFieldTemplates.C"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
static wordHashSet supportedTypes()
{
    // typedef DimensionedField<Type, areaMesh> DimFieldType;
    typedef GeometricField<Type, faPatchField, areaMesh> AreaFieldType;

    wordHashSet types;

    // TODO: types.insert(DimFieldType::typeName);
    types.insert(AreaFieldType::typeName);

    return types;
}

} // end of Foam


Foam::wordHashSet Foam::vtk::faMeshAdaptor::knownFields
(
    const wordRes& selectFields
) const
{
    wordHashSet allFields;

    // Quick exit if no volume fields can be converted.
    // This could be refined
    HashTable<wordHashSet> objects = mesh_.mesh().classes(selectFields);

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


void Foam::vtk::faMeshAdaptor::convertAreaFields
(
    const wordRes& selectFields
)
{
    // Quick exit if no fields can be converted.
    // This could be refined
    HashTable<wordHashSet> objects = mesh_.mesh().classes(selectFields);

    if (objects.empty())
    {
        return;
    }

    convertAreaFields<scalar>(selectFields);
    convertAreaFields<vector>(selectFields);
    convertAreaFields<sphericalTensor>(selectFields);
    convertAreaFields<symmTensor>(selectFields);
    convertAreaFields<tensor>(selectFields);
}


// ************************************************************************* //
