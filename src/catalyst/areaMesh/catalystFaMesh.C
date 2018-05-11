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

#include "catalystFaMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "faMesh.H"
#include "fvMesh.H"

#include <vtkCPDataDescription.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace catalyst
{
    defineTypeNameAndDebug(faMeshInput, 0);
    addNamedToRunTimeSelectionTable
    (
        catalystInput,
        faMeshInput,
        dictionary,
        area
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::catalyst::faMeshInput::update()
{
    // Backend requires a corresponding mesh
    backends_.filterKeys
    (
        [this](const word& k){ return meshes_.found(k); }
    );

    forAllConstIters(meshes_, iter)
    {
        if (!backends_.found(iter.key()))
        {
            backends_.set
            (
                iter.key(),
                new Foam::vtk::faMeshAdaptor(*(iter.object()))
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::catalyst::faMeshInput::faMeshInput
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    catalystInput(name),
    time_(runTime),
    regionName_(),
    selectAreas_(),
    selectFields_(),
    meshes_(),
    backends_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::catalyst::faMeshInput::read(const dictionary& dict)
{
    catalystInput::read(dict);

    selectAreas_.clear();
    selectFields_.clear();
    backends_.clear();

    regionName_ =
        dict.lookupOrDefault<word>("region", polyMesh::defaultRegion);

    const objectRegistry& obr =
        time_.lookupObject<objectRegistry>(regionName_);

    // All possible meshes for the given region
    meshes_ = obr.lookupClass<faMesh>();

    dict.readIfPresent("areas", selectAreas_);

    if (selectAreas_.empty())
    {
        word areaName;
        if (!dict.readIfPresent("area", areaName))
        {
            wordList available = obr.sortedNames<faMesh>();

            if (available.size())
            {
                areaName = available.first();
            }
        }

        if (!areaName.empty())
        {
            selectAreas_.resize(1);
            selectAreas_.first() = areaName;
        }
    }

    // Restrict to specified meshes
    meshes_.filterKeys(selectAreas_);

    dict.lookup("fields") >> selectFields_;

    return true;
}


void Foam::catalyst::faMeshInput::update(polyMesh::readUpdateState state)
{
    // Trigger change of state

    const objectRegistry& obr =
        time_.lookupObject<objectRegistry>(regionName_);

    // Be really paranoid and verify if the mesh actually exists
    const wordList regionNames(backends_.toc());

    for (const word& regionName : regionNames)
    {
        if (meshes_.found(regionName) && obr.found(regionName))
        {
            backends_[regionName]->updateState(state);
        }
        else
        {
            backends_.erase(regionName);
            meshes_.erase(regionName);
        }
    }
}


Foam::label Foam::catalyst::faMeshInput::addChannels(dataQuery& dataq)
{
    update();   // Enforce sanity for backends and adaptor

    if (backends_.empty())
    {
        return 0;
    }

    // Gather all fields that we know how to convert
    wordHashSet allFields;
    forAllConstIters(backends_, iter)
    {
        allFields += iter.object()->knownFields(selectFields_);
    }


    dataq.set(name(), allFields);

    return 1;
}


bool Foam::catalyst::faMeshInput::convert
(
    dataQuery& dataq,
    outputChannels& outputs
)
{
    const wordList regionNames(backends_.sortedToc());

    if (regionNames.empty())
    {
        return false;  // skip - not available
    }

    // Single channel only

    label nChannels = 0;

    if (dataq.found(name()))
    {
        ++nChannels;
    }

    if (!nChannels)
    {
        return false;  // skip - not requested
    }


    // TODO: currently don't rely on the results from expecting much at all

    // Each region in a separate block
    unsigned int blockNo = 0;
    for (const word& regionName : regionNames)
    {
        auto dataset =
            backends_[regionName]->output(selectFields_);

        {
            const fileName channel = name();

            if (dataq.found(channel))
            {
                // Get existing or new
                vtkSmartPointer<vtkMultiBlockDataSet> block =
                    outputs.lookup
                    (
                        channel,
                        vtkSmartPointer<vtkMultiBlockDataSet>::New()
                    );

                block->SetBlock(blockNo, dataset);

                block->GetMetaData(blockNo)->Set
                (
                    vtkCompositeDataSet::NAME(),
                    regionName
                );

                outputs.set(channel, block);  // overwrites existing
            }
        }

        ++blockNo;
    }

    return true;
}


Foam::Ostream& Foam::catalyst::faMeshInput::print(Ostream& os) const
{
    os  << name() << nl
        <<"    areas   " << flatOutput(selectAreas_) << nl
        <<"    meshes  " << flatOutput(meshes_.sortedToc()) << nl
        <<"    fields  " << flatOutput(selectFields_) << nl;

    return os;
}


// ************************************************************************* //
