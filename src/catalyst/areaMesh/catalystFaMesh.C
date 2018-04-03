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
#include "catalystCoprocess.H"
#include "addToRunTimeSelectionTable.H"
#include "faMesh.H"
#include "fvMesh.H"

#include <vtkNew.h>
#include <vtkCPDataDescription.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkInformation.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(catalystFaMesh, 0);
    addToRunTimeSelectionTable(functionObject, catalystFaMesh, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::catalystFaMesh::updateState
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

Foam::functionObjects::catalystFaMesh::catalystFaMesh
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    selectAreas_(),
    selectFields_(),
    scripts_(),
    meshes_(),
    backends_(),
    adaptor_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::catalystFaMesh::~catalystFaMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::catalystFaMesh::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    int debugLevel = 0;
    if (dict.readIfPresent("debug", debugLevel))
    {
        catalystCoprocess::debug = debugLevel;
    }

    // All possible meshes
    meshes_ = mesh_.lookupClass<faMesh>();

    selectAreas_.clear();
    dict.readIfPresent("areas", selectAreas_);

    if (selectAreas_.empty())
    {
        word areaName;

        if (!dict.readIfPresent("area", areaName))
        {
            wordList available = mesh_.sortedNames<faMesh>();

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
    meshes_.filterKeys(wordRes(selectAreas_));

    dict.lookup("fields") >> selectFields_;
    dict.lookup("scripts") >> scripts_;         // Python scripts

    catalystCoprocess::expand(scripts_, dict);  // Expand and check availability

    Info<< type() << " " << name() << ":" << nl
        <<"    areas   " << flatOutput(selectAreas_) << nl
        <<"    meshes  " << flatOutput(meshes_.sortedToc()) << nl
        <<"    fields  " << flatOutput(selectFields_) << nl
        <<"    scripts " << scripts_ << nl;

    // Run-time modification of pipeline
    if (adaptor_.valid())
    {
        adaptor_().reset(scripts_);
    }

    // Ensure consistency - only retain backends with corresponding mesh region
    backends_.retain(meshes_);

    return true;
}


bool Foam::functionObjects::catalystFaMesh::execute()
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
                    new Foam::vtk::faMeshAdaptor(*(iter.object()))
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
        vtk::faMeshAdaptor::channelNames.names(),
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

    auto output = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    // TODO: currently don't rely on the results from expecting much at all

    // Each region in a separate block
    unsigned int regionNo = 0;
    for (const word& regionName : regionNames)
    {
        auto pieces = backends_[regionName]->output(selectFields_);

        output->SetBlock(regionNo, pieces);

        output->GetMetaData(regionNo)->Set
        (
            vtkCompositeDataSet::NAME(),
            regionName
        );
        ++regionNo;
    }

    if (regionNo)
    {
        Log << type() << ": send data" << nl;

        adaptor_().process(dataq, output);
    }

    return true;
}


bool Foam::functionObjects::catalystFaMesh::write()
{
    return true;
}


void Foam::functionObjects::catalystFaMesh::updateMesh(const mapPolyMesh&)
{
    updateState(polyMesh::TOPO_CHANGE);
}


void Foam::functionObjects::catalystFaMesh::movePoints(const polyMesh&)
{
    updateState(polyMesh::POINTS_MOVED);
}


// ************************************************************************* //
