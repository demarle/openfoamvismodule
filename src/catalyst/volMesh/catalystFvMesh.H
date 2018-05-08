/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
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

Class
    Foam::functionObjects::catalystFvMesh

Group
    grpUtilitiesFunctionObjects

Description
    A Paraview Catalyst adaptor for OpenFOAM fvMesh regions.

    The output comprises up to three channels ("input", "mesh", "patches"),
    each of which is a multi-block dataset.

    Example of function object specification:
    \verbatim
    catalyst
    {
        type            catalyst;
        libs            ("libcatalystFoam.so");
        regions         ( ".*Solid" )
        fields          (U p);
        scripts         ( ... );
        executeControl  timeStep;
        executeInterval 1;
    }
    \endverbatim

Usage
    \table
        Property     | Description                 | Required    | Default
        type         | catalyst                    | yes         |
        log          | report extra information    | no          | false
        mkdir        | initial directory to create | no          |
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
    Foam::vtk::fvMeshAdaptor

SourceFiles
    catalystFvMesh.C
    catalystFvMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_catalystFvMesh_H
#define functionObjects_catalystFvMesh_H

#include "functionObject.H"
#include "foamVtkFvMeshAdaptor.H"
#include "wordRes.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class catalystCoprocess;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class catalystFvMesh Declaration
\*---------------------------------------------------------------------------*/

class catalystFvMesh
:
    public functionObject
{
protected:

    // Protected data

        //- Reference to the time database
        const Time& time_;

        //- The output directory
        fileName outputDir_;

        //- Python scripts for the catalyst pipeline
        stringList scripts_;

        //- The catalyst coprocess
        autoPtr<catalystCoprocess> adaptor_;


        //- Requested names of regions to process
        wordRes selectRegions_;

        //- Names of fields to process
        wordRes selectFields_;

        //- Pointers to the requested mesh regions
        HashTable<const fvMesh*> meshes_;

        //- Backends for OpenFOAM to VTK translation (with internal caching)
        HashPtrTable<vtk::fvMeshAdaptor> backends_;


    // Protected Member Functions

        //- Common boilerplate settings
        bool readBasics(const dictionary& dict);

        //- On movement
        void updateState(polyMesh::readUpdateState state);


        //- No copy construct
        catalystFvMesh(const catalystFvMesh&) = delete;

        //- No copy assignment
        void operator=(const catalystFvMesh&) = delete;

public:

    //- Runtime type information
    TypeName("catalyst");


    // Constructors

        //- Construct from Time and dictionary
        catalystFvMesh
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~catalystFvMesh();


    // Member Functions

        //- Read the specification
        virtual bool read(const dictionary& dict);

        //- Execute catalyst pipelines
        virtual bool execute();

        //- Write - does nothing
        virtual bool write();

        //- On end - provide feedback about disconnecting from catatyst.
        virtual bool end();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for mesh point-motion
        virtual void movePoints(const polyMesh& mesh);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
