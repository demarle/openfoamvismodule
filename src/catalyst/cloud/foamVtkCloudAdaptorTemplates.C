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

#ifndef foamVtkCloudAdaptorTemplates_C
#define foamVtkCloudAdaptorTemplates_C

#include "foamVtkTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class UnaryMatchPredicate>
Foam::label Foam::vtk::cloudAdaptor::convertLagrangianFields
(
    vtkPolyData* vtkmesh,
    const objectRegistry& obr,
    const UnaryMatchPredicate& pred
)
{
    label nFields = 0;

    wordHashSet names(obr.names<IOField<Type>>());

    // Restrict to specified names
    names.filterKeys(pred);

    // Avoid converting particle positions as a field too.
    names.filterKeys
    (
        [](const word& k)
        {
            return k.startsWith("position") || k.startsWith("coordinate");
        },
        true  // prune
    );

    const wordList fieldNames(names.sortedToc());

    for (const word& fieldName : fieldNames)
    {
        const auto& fld = obr.lookupObject<IOField<Type>>(fieldName);

        vtkSmartPointer<vtkFloatArray> data =
            vtk::Tools::convertFieldToVTK
            (
                fieldName,
                fld
            );

        // Provide identical data as cell and as point data
        vtkmesh->GetCellData()->AddArray(data);
        vtkmesh->GetPointData()->AddArray(data);

        ++nFields;
    }

    return nFields;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
