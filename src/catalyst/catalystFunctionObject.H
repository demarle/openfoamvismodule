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

Class
    Foam::functionObjects::catalystFunctionObject

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
        scripts         ( ... );
        executeControl  timeStep;
        executeInterval 1;

        inputs
        {
            type        default;
            regions     ( ".*Solid" )
            fields      (U p);
        }
    }
    \endverbatim

Usage
    \table
        Property     | Description                 | Required    | Default
        type         | catalyst                    | yes         |
        log          | report extra information    | no          | false
        mkdir        | initial directory to create | no          |
        scripts      | Python pipeline scripts     | yes         |
    \endtable

Note
    The execution frequency can be defined by the functionObject and
    by the Catalyst pipeline.

See also
    Foam::functionObjects::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::timeControl
    Foam::catalyst::coprocess
    Foam::catalyst::cloudInput
    Foam::catalyst::faMeshInput
    Foam::catalyst::fvMeshInput

SourceFiles
    catalystFunctionObject.C
    catalystFunctionObjectTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_catalystFunctionObject_H
#define functionObjects_catalystFunctionObject_H

#include "className.H"
#include "wordList.H"
#include "stringList.H"
#include "polyMesh.H"
#include "PtrList.H"
#include "functionObject.H"
#include "catalystCoprocess.H"
#include "catalystInput.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class catalystFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class catalystFunctionObject
:
    public functionObject
{
    // Private data

        //- Reference to the time database
        const Time& time_;

        //- The output directory
        fileName outputDir_;

        //- Python scripts for the catalyst pipeline
        stringList scripts_;

        //- The catalyst coprocess
        autoPtr<catalyst::coprocess> adaptor_;

        //- Pointers to the requested mesh regions
        PtrList<catalyst::catalystInput> inputs_;


    // Private Member Functions

        //- No copy construct
        catalystFunctionObject(const catalystFunctionObject&) = delete;

        //- No copy assignment
        void operator=(const catalystFunctionObject&) = delete;

public:

    // Static Methods

        //- Expand strings as filenames, retaining only those that exist
        static label expand(List<string>& scripts, const dictionary& dict);


    //- Runtime type information
    TypeName("catalyst");


    // Constructors

        //- Construct from Time and dictionary
        catalystFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~catalystFunctionObject();


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