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

Class
    Foam::functionObjects::catalystCloud

Group
    grpUtilitiesFunctionObjects

Description
    A Paraview Catalyst adaptor for a single OpenFOAM cloud (lagrangian).

    The output comprises a single "cloud" channel, which is a multi-block
    dataset (one block per cloud).

    Example of function object specification:
    \verbatim
    catalystCloud1
    {
        type            catalyst::cloud;
        libs            ("libcatalystFoam.so");
        cloud           NAME;
        fields          (U T rho);
        scripts         ( ... );
        executeControl  timeStep;
        executeInterval 1;
    }
    \endverbatim

Usage
    \table
        Property     | Description                 | Required    | Default value
        type         | catalyst::cloud             | yes         |
        log          | report extra information    | no          | false
        mkdir        | initial directory to create | no          |
        cloud        |                             | no          | defaultCloud
        clouds       | wordRe list of clouds       | no          |
        region       |                             | no          | region0
        regions      | wordRe list of regions      | no          |
        fields       | wordRe list of fields       | yes         |
        scripts      | Python pipeline scripts     | yes         |
    \endtable

Note
    The execution frequency can be defined by the functionObject and
    by the Catalyst pipeline.

See also
    Foam::functionObjects::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl
    Foam::catalystCoprocess
    Foam::vtk::cloudAdaptor

SourceFiles
    catalystCloud.C
    catalystCloudTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_catalystCloud_H
#define functionObjects_catalystCloud_H

#include "fvMeshFunctionObject.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class catalystCoprocess;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class catalystCloud Declaration
\*---------------------------------------------------------------------------*/

class catalystCloud
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- Requested names of clouds to process
        wordRes selectClouds_;

        //- Subset of cloud fields to process
        wordRes selectFields_;

        //- Python scripts for the catalyst pipeline
        stringList scripts_;

        //- The catalyst coprocess
        autoPtr<catalystCoprocess> adaptor_;


    // Protected Member Functions

        //- Common boilerplate settings
        bool readBasics(const dictionary& dict);


        //- No copy construct
        catalystCloud(const catalystCloud&) = delete;

        //- No copy assignment
        void operator=(const catalystCloud&) = delete;

public:

    //- Runtime type information
    TypeName("catalyst::cloud");


    // Constructors

        //- Construct from Time and dictionary
        catalystCloud
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~catalystCloud();


    // Member Functions

        //- Read the specification
        virtual bool read(const dictionary& dict);

        //- Execute catalyst pipelines
        virtual bool execute();

        //- Write - does nothing
        virtual bool write();

        //- On end - provide feedback about disconnecting from catatyst.
        virtual bool end();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
