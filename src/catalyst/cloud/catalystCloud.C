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

#include "catalystCloud.H"
#include "catalystCoprocess.H"
#include "cloud.H"
#include "foamVtkCloudAdaptor.H"
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
    defineTypeNameAndDebug(catalystCloud, 0);
    addToRunTimeSelectionTable(functionObject, catalystCloud, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::catalystCloud::readBasics(const dictionary& dict)
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::catalystCloud::catalystCloud
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    selectClouds_(),
    selectFields_(),
    adaptor_()
{
    if (postProcess)
    {
        // Disable for post-process mode.
        // Emit as FatalError for the try/catch in the caller.
        FatalError
            << type() << " disabled in post-process mode"
            << exit(FatalError);
    }
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::catalystCloud::~catalystCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::catalystCloud::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    readBasics(dict);

    selectClouds_.clear();
    dict.readIfPresent("clouds", selectClouds_);

    if (selectClouds_.empty())
    {
        selectClouds_.resize(1);
        selectClouds_.first() =
            dict.lookupOrDefault<word>("cloud", cloud::defaultName);
    }

    selectFields_.clear();
    dict.readIfPresent("fields", selectFields_);


    Info<< type() << " " << name() << ":" << nl
        <<"    clouds  " << flatOutput(selectClouds_) << nl
        <<"    fields  " << flatOutput(selectFields_) << nl
        <<"    scripts " << scripts_ << nl;

    if (adaptor_.valid())
    {
        // Run-time modification of pipeline
        adaptor_().reset(scripts_);
    }

    return true;
}


bool Foam::functionObjects::catalystCloud::execute()
{
    const wordList cloudNames(mesh_.sortedNames<cloud>(selectClouds_));

    if (cloudNames.empty())
    {
        return true;
    }

    // Enforce sanity for backends and adaptor
    {
        if (!adaptor_.valid())
        {
            adaptor_.reset(new catalystCoprocess());
            adaptor_().reset(scripts_);
        }
    }


    // Difficult to get the names of the fields from particles
    // ... need to skip for now
    wordHashSet allFields;


    // Data description for co-processing
    vtkNew<vtkCPDataDescription> descrip;

    // Form data query for catalyst
    catalystCoprocess::dataQuery dataq
    (
        vtk::cloudAdaptor::channelNames.names(),
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

    // Each cloud in a separate block.
    unsigned int cloudNo = 0;
    for (const word& cloudName : cloudNames)
    {
        auto pieces =
            vtk::cloudAdaptor(mesh_).getCloud(cloudName, selectFields_);

        output->SetBlock(cloudNo, pieces);
        output->GetMetaData(cloudNo)->Set
        (
            vtkCompositeDataSet::NAME(),
            cloudName
        );

        ++cloudNo;
    }

    if (cloudNo)
    {
        Log << type() << ": send data" << nl;

        adaptor_().process(dataq, output);
    }

    return true;
}


bool Foam::functionObjects::catalystCloud::write()
{
    return true;
}


bool Foam::functionObjects::catalystCloud::end()
{
    // Only here for extra feedback
    if (log && adaptor_.valid())
    {
        Info<< type() << ": Disconnecting ParaView Catalyst..." << nl;
    }

    adaptor_.clear();
    return true;
}


// ************************************************************************* //
