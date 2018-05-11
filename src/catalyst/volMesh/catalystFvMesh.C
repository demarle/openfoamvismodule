/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2018 CINECA
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

#include "catalystFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "fileNameList.H"

#include <vtkCPDataDescription.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace catalyst
{
    defineTypeNameAndDebug(fvMeshInput, 0);
    addNamedToRunTimeSelectionTable
    (
        catalystInput,
        fvMeshInput,
        dictionary,
        default
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::catalyst::fvMeshInput::channelName(channelEnum chan) const
{
    if (chan == channelEnum::INPUT)
    {
        return name();
    }

    return name()/vtk::fvMeshAdaptor::channelNames[chan];
}


void Foam::catalyst::fvMeshInput::update()
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
                new Foam::vtk::fvMeshAdaptor(*(iter.object()), decompose_)
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::catalyst::fvMeshInput::fvMeshInput
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    catalystInput(name),
    time_(runTime),
    selectRegions_(),
    selectFields_(),
    decompose_(false),
    meshes_(),
    backends_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::catalyst::fvMeshInput::read(const dictionary& dict)
{
    catalystInput::read(dict);

    backends_.clear();
    selectFields_.clear();
    selectRegions_.clear();
    decompose_ = dict.lookupOrDefault("decompose", false);

    // All possible meshes
    meshes_ = time_.lookupClass<fvMesh>();

    dict.readIfPresent("regions", selectRegions_);

    if (selectRegions_.empty())
    {
        selectRegions_.resize(1);
        selectRegions_.first() =
            dict.lookupOrDefault<word>("region", polyMesh::defaultRegion);
    }

    // Restrict to specified meshes
    meshes_.filterKeys(selectRegions_);

    dict.lookup("fields") >> selectFields_;

    return true;
}


void Foam::catalyst::fvMeshInput::update(polyMesh::readUpdateState state)
{
    // Trigger change of state

    // Be really paranoid and verify if the mesh actually exists
    const wordList regionNames(backends_.toc());

    for (const word& regionName : regionNames)
    {
        if (meshes_.found(regionName) && time_.found(regionName))
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


Foam::label Foam::catalyst::fvMeshInput::addChannels(dataQuery& dataq)
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


    // This solution may be a temporary measure...

    dataq.set(name(), allFields);  // channel::INPUT
    dataq.set(channelName(channelEnum::MESH),    allFields);
    dataq.set(channelName(channelEnum::PATCHES), allFields);

    return 3;
}


bool Foam::catalyst::fvMeshInput::convert
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

    // Multi-channel

    // This solution may be a temporary measure...

    unsigned whichChannels = channelEnum::NONE;

    // channel::INPUT
    if (dataq.found(name()))
    {
        whichChannels |= channelEnum::INPUT;
    }

    // channel::MESH
    if (dataq.found(channelName(channelEnum::MESH)))
    {
        whichChannels |= channelEnum::MESH;
    }

    // channel::PATCHES
    if (dataq.found(channelName(channelEnum::PATCHES)))
    {
        whichChannels |= channelEnum::PATCHES;
    }

    if (channelEnum::NONE == whichChannels)
    {
        return false;  // skip - not requested
    }


    // TODO: currently don't rely on the results from expecting much at all

    // Each region goes into a separate block.
    unsigned int blockNo = 0;
    for (const word& regionName : regionNames)
    {
        // define/redefine output channels
        backends_[regionName]->channels(whichChannels);

        vtkSmartPointer<vtkMultiBlockDataSet> dataset =
            backends_[regionName]->output(selectFields_);

        // MESH = block 0
        {
            const unsigned int channelNo = 0;

            const fileName channel(channelName(channelEnum::MESH));

            if (dataq.found(channel))
            {
                // Get existing or new
                vtkSmartPointer<vtkMultiBlockDataSet> block =
                    outputs.lookup
                    (
                        channel,
                        vtkSmartPointer<vtkMultiBlockDataSet>::New()
                    );

                block->SetBlock(blockNo, dataset->GetBlock(channelNo));

                block->GetMetaData(blockNo)->Set
                (
                    vtkCompositeDataSet::NAME(),
                    regionName
                );

                outputs.set(channel, block); // overwrites existing
            }
        }

        // PATCHES = block 1
        {
            const unsigned int channelNo = 1;

            const fileName channel(channelName(channelEnum::PATCHES));

            if (dataq.found(channel))
            {
                // Get existing or new
                vtkSmartPointer<vtkMultiBlockDataSet> block =
                    outputs.lookup
                    (
                        channel,
                        vtkSmartPointer<vtkMultiBlockDataSet>::New()
                    );

                block->SetBlock(blockNo, dataset->GetBlock(channelNo));

                block->GetMetaData(blockNo)->Set
                (
                    vtkCompositeDataSet::NAME(),
                    regionName
                );

                outputs.set(channel, block);  // overwrites existing
            }
        }

        // INPUT
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


Foam::Ostream& Foam::catalyst::fvMeshInput::print(Ostream& os) const
{
    os  << name() << nl
        <<"    regions " << flatOutput(selectRegions_) << nl
        <<"    meshes  " << flatOutput(meshes_.sortedToc()) << nl
        <<"    fields  " << flatOutput(selectFields_) << nl;

    return os;
}


// ************************************************************************* //
