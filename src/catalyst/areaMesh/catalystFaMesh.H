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
    Foam::functionObjects::catalystFaMesh

Group
    grpUtilitiesFunctionObjects

Description
    A Paraview Catalyst adaptor for OpenFOAM faMesh regions.

    The output comprises a single "input" channel, which is a multi-block
    dataset (one block per area mesh).

    Example of function object specification:
    \verbatim
    catalyst
    {
        type            catalyst::area;
        libs            ("libcatalystFoam.so");
        fields          (U p);
        scripts         ( ... );
        executeControl  timeStep;
        executeInterval 1;
    }
    \endverbatim

Usage
    \table
        Property     | Description                 | Required    | Default
        type         | catalyst::area              | yes         |
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
    Foam::vtk::faMeshAdaptor

SourceFiles
    catalystFaMesh.C
    catalystFaMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_catalystFaMesh_H
#define functionObjects_catalystFaMesh_H

#include "fvMeshFunctionObject.H"
#include "foamVtkFaMeshAdaptor.H"
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
                          Class catalystFaMesh Declaration
\*---------------------------------------------------------------------------*/

class catalystFaMesh
:
    public fvMeshFunctionObject
{
protected:

    // Protected data

        //- The output directory
        fileName outputDir_;

        //- Python scripts for the catalyst pipeline
        stringList scripts_;

        //- The catalyst coprocess
        autoPtr<catalystCoprocess> adaptor_;


        //- Requested names of areas to process
        wordRes selectAreas_;

        //- Names of fields to process
        wordRes selectFields_;

        //- Pointers to the requested mesh regions
        HashTable<const faMesh*> meshes_;

        //- Backends for OpenFOAM to VTK translation (with internal caching)
        HashPtrTable<vtk::faMeshAdaptor> backends_;


    // Protected Member Functions

        //- Common boilerplate settings
        bool readBasics(const dictionary& dict);

        //- On movement
        void updateState(polyMesh::readUpdateState state);


        //- No copy construct
        catalystFaMesh(const catalystFaMesh&) = delete;

        //- No copy assignment
        void operator=(const catalystFaMesh&) = delete;

public:

    //- Runtime type information
    TypeName("catalyst::area");


    // Constructors

        //- Construct from Time and dictionary
        catalystFaMesh
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~catalystFaMesh();


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
