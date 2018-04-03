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

#ifndef foamVtkFaMeshAdaptorFieldTemplates_C
#define foamVtkFaMeshAdaptorFieldTemplates_C

// OpenFOAM includes
#include "error.H"
#include "emptyFvPatchField.H"
#include "wallPolyPatch.H"
#include "volPointInterpolation.H"
#include "zeroGradientFvPatchField.H"
#include "interpolatePointToCell.H"

// vtk includes
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::faMeshAdaptor::convertAreaFields
(
    const wordRes& selectFields
)
{
    typedef GeometricField<Type, faPatchField, areaMesh> FieldType;

    // Restrict to GeometricField<Type, ...>
    const wordList names(mesh_.mesh().sortedNames<FieldType>(selectFields));

    for (const word& fieldName : names)
    {
        convertAreaField
        (
            mesh_.mesh().lookupObject<FieldType>(fieldName)
        );
    }
}


template<class Type>
void Foam::vtk::faMeshAdaptor::convertAreaField
(
    const GeometricField<Type, faPatchField, areaMesh>& fld
)
{
    const auto& longName = internalName;

    auto iter = cachedVtp_.find(longName);
    if (!iter.found() || !iter.object().dataset)
    {
        // Should not happen, but for safety require a vtk geometry
        return;
    }
    foamVtpData& vtpData = iter.object();
    auto dataset = vtpData.dataset;

    vtkSmartPointer<vtkFloatArray> cdata = convertAreaFieldToVTK
    (
        fld,
        vtpData
    );
    dataset->GetCellData()->AddArray(cdata);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// low-level conversions
//

template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::faMeshAdaptor::convertAreaFieldToVTK
(
    const GeometricField<Type, faPatchField, areaMesh>& fld,
    const foamVtpData& vtpData
) const
{
    const int nComp(pTraits<Type>::nComponents);

    auto data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetName(fld.name().c_str());
    data->SetNumberOfComponents(nComp);
    data->SetNumberOfTuples(fld.size());

    const label len = fld.size();

    if (debug)
    {
        Info<< "convert areaField: "
            << fld.name()
            << " size=" << len
            << " nComp=" << nComp << endl;
    }

    float scratch[nComp];
    for (label i=0; i < len; ++i)
    {
        const Type& t = fld[i];
        for (direction d=0; d<nComp; ++d)
        {
            scratch[d] = component(t, d);
        }
        remapTuple<Type>(scratch);

        data->SetTuple(i, scratch);
    }

    return data;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
