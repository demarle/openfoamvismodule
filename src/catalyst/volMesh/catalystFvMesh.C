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
#include "catalystCoprocess.H"
#include "addToRunTimeSelectionTable.H"

#include <vtkNew.h>
#include <vtkCPDataDescription.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(catalystFvMesh, 0);
    addToRunTimeSelectionTable(functionObject, catalystFvMesh, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::catalystFvMesh::readBasics(const dictionary& dict)
{
    int debugLevel = 0;
    if (dict.readIfPresent("debug", debugLevel))
    {
        catalystCoprocess::debug = debugLevel;
    }

    fileName outputDir;
    if (dict.readIfPresent("mkdir", outputDir))
    {
        outputDir.expand();
        outputDir.clean();
        Foam::mkDir(outputDir);
    }

    dict.lookup("scripts") >> scripts_;         // Python scripts
    catalystCoprocess::expand(scripts_, dict);  // Expand and check availability

    return true;
}


void Foam::functionObjects::catalystFvMesh::updateState
(
    polyMesh::readUpdateState state
)
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::catalystFvMesh::catalystFvMesh
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    selectRegions_(),
    selectFields_(),
    scripts_(),
    meshes_(),
    backends_(),
    adaptor_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::catalystFvMesh::~catalystFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::catalystFvMesh::read(const dictionary& dict)
{
    functionObject::read(dict);

    // Common settings
    int debugLevel = 0;
    if (dict.readIfPresent("debug", debugLevel))
    {
        catalystCoprocess::debug = debugLevel;
    }

    fileName outputDir;
    if (dict.readIfPresent("mkdir", outputDir))
    {
        outputDir.expand();
        outputDir.clean();
        Foam::mkDir(outputDir);
    }

    dict.lookup("scripts") >> scripts_;         // Python scripts
    catalystCoprocess::expand(scripts_, dict);  // Expand and check availability


    // All possible meshes
    meshes_ = time_.lookupClass<fvMesh>();

    selectRegions_.clear();
    dict.readIfPresent("regions", selectRegions_);

    if (selectRegions_.empty())
    {
        selectRegions_.resize(1);
        selectRegions_.first() =
            dict.lookupOrDefault<word>("region", polyMesh::defaultRegion);
    }

    // Restrict to specified meshes
    meshes_.filterKeys(wordRes(selectRegions_));

    dict.lookup("fields") >> selectFields_;


    Info<< type() << " " << name() << ":" << nl
        <<"    regions " << flatOutput(selectRegions_) << nl
        <<"    meshes  " << flatOutput(meshes_.sortedToc()) << nl
        <<"    fields  " << flatOutput(selectFields_) << nl
        <<"    scripts " << scripts_ << nl;

    if (adaptor_.valid())
    {
        // Run-time modification of pipeline
        adaptor_().reset(scripts_);
    }

    // Ensure consistency - only retain backends with corresponding mesh region
    backends_.retain(meshes_);

    return true;
}


bool Foam::functionObjects::catalystFvMesh::execute()
{
    const wordList regionNames(meshes_.sortedToc());

    if (regionNames.empty())
    {
        return false;
    }

    // Enforce sanity for backends and adaptor
    {
        bool updateAdaptor = false;
        forAllConstIters(meshes_, iter)
        {
            if (!backends_.found(iter.key()))
            {
                backends_.insert
                (
                    iter.key(),
                    new Foam::vtk::fvMeshAdaptor(*(iter.object()))
                );
                updateAdaptor = true;
            }
        }

        if (updateAdaptor && !adaptor_.valid())
        {
            adaptor_.reset(new catalystCoprocess());
            adaptor_().reset(scripts_);
        }
    }


    // Gather all fields that we know how to convert

    wordHashSet allFields;
    forAllConstIters(backends_, iter)
    {
        allFields += iter.object()->knownFields(selectFields_);
    }


    // Data description for co-processing
    vtkNew<vtkCPDataDescription> descrip;

    // Form data query for catalyst
    catalystCoprocess::dataQuery dataq
    (
        vtk::fvMeshAdaptor::channelNames.names(),
        time_,  // timeQuery
        descrip.Get()
    );

    // Query catalyst
    const HashTable<wordHashSet> expecting(adaptor_().query(dataq, allFields));

    if (catalystCoprocess::debug)
    {
        if (expecting.empty())
        {
            Info<< type() << ": expecting no data" << nl;
        }
        else
        {
            Info<< type() << ": expecting data " << expecting << nl;
        }
    }

    if (expecting.empty())
    {
        return true;
    }

    HashTable<vtkSmartPointer<vtkMultiBlockDataSet>> outputs;

    // TODO: currently don't rely on the results from expecting much at all

    // Each region in a separate block.
    unsigned int regionNo = 0;
    for (const word& regionName : regionNames)
    {
        // (re)define output channels
        backends_[regionName]->channels(expecting.toc());

        vtkSmartPointer<vtkMultiBlockDataSet> dataset =
            backends_[regionName]->output(selectFields_);

        {
            const unsigned int channelNo = 0; // MESH

            const word& channelName =
                vtk::fvMeshAdaptor::channelNames
                [vtk::fvMeshAdaptor::channel::MESH];

            if (expecting.found(channelName))
            {
                // Get existing or new
                vtkSmartPointer<vtkMultiBlockDataSet> block =
                    outputs.lookup
                    (
                        channelName,
                        vtkSmartPointer<vtkMultiBlockDataSet>::New()
                    );

                block->SetBlock(regionNo, dataset->GetBlock(channelNo));

                block->GetMetaData(regionNo)->Set
                (
                    vtkCompositeDataSet::NAME(),
                    regionName
                );

                outputs.set(channelName, block);  // overwrite existing
            }
        }

        {
            const unsigned int channelNo = 1; // PATCHES

            const word& channelName =
                vtk::fvMeshAdaptor::channelNames
                [vtk::fvMeshAdaptor::channel::PATCHES];

            if (expecting.found(channelName))
            {
                // Get existing or new
                vtkSmartPointer<vtkMultiBlockDataSet> block =
                    outputs.lookup
                    (
                        channelName,
                        vtkSmartPointer<vtkMultiBlockDataSet>::New()
                    );

                block->SetBlock(regionNo, dataset->GetBlock(channelNo));

                block->GetMetaData(regionNo)->Set
                (
                    vtkCompositeDataSet::NAME(),
                    regionName
                );

                outputs.set(channelName, block);  // overwrite existing
            }
        }

        {
            const word& channelName =
                vtk::fvMeshAdaptor::channelNames
                [vtk::fvMeshAdaptor::channel::INPUT];

            if (expecting.found(channelName))
            {
                // Get existing or new
                vtkSmartPointer<vtkMultiBlockDataSet> block =
                    outputs.lookup
                    (
                        channelName,
                        vtkSmartPointer<vtkMultiBlockDataSet>::New()
                    );

                block->SetBlock(regionNo, dataset);

                block->GetMetaData(regionNo)->Set
                (
                    vtkCompositeDataSet::NAME(),
                    regionName
                );

                outputs.set(channelName, block);  // overwrite existing
            }
        }

        ++regionNo;
    }

    if (regionNo)
    {
        Log << type() << ": send data" << nl;

        adaptor_().process(dataq, outputs);
    }

    return true;
}


bool Foam::functionObjects::catalystFvMesh::write()
{
    return true;
}


bool Foam::functionObjects::catalystFvMesh::end()
{
    // Only here for extra feedback
    if (log && adaptor_.valid())
    {
        Info<< type() << ": Disconnecting ParaView Catalyst..." << nl;
    }

    adaptor_.clear();
    return true;
}


void Foam::functionObjects::catalystFvMesh::updateMesh(const mapPolyMesh&)
{
    updateState(polyMesh::TOPO_CHANGE);
}


void Foam::functionObjects::catalystFvMesh::movePoints(const polyMesh&)
{
    updateState(polyMesh::POINTS_MOVED);
}


// ************************************************************************* //
