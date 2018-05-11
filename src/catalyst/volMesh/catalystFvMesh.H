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
    Foam::catalyst::fvMeshInput

Description
    A Paraview Catalyst source for OpenFOAM fvMesh regions.

    The source comprises up to three internal channels
    ("input", "mesh", "patches"), each of which is a multi-block dataset.
    On output, the "input" sub-channel receives the name of the source,
    while the "mesh" and "patches" sub-channels are prefixed with the
    name of the source.

    Example specification:
    \verbatim
    myName
    {
        type            default;
        regions         ( ".*Solid" )
        fields          (U p);
    }
    \endverbatim

    Resulting in three output channels:
    "myName", "myName/mesh", "myName/patches"

Usage
    \table
        Property     | Description                      | Required    | Default
        type         | default                          | no          | default
        region       |                                  | no          | region0
        regions      | wordRe list of regions           | no          |
        fields       | wordRe list of fields            | yes         |
        decompose    | decompose polyhedra              | no          | false
    \endtable

See also
    Foam::catalyst::catalystInput
    Foam::vtk::fvMeshAdaptor

SourceFiles
    catalystFvMesh.C
    catalystFvMeshTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef catalyst_fvMeshInput_H
#define catalyst_fvMeshInput_H

#include "wordRes.H"
#include "HashPtrTable.H"
#include "catalystInput.H"
#include "foamVtkFvMeshAdaptor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace catalyst
{

/*---------------------------------------------------------------------------*\
                    Class catalyst::fvMeshInput Declaration
\*---------------------------------------------------------------------------*/

class fvMeshInput
:
    public catalystInput
{
protected:

    // Protected data

        //- Reference to the time database
        const Time& time_;

        //- Requested names of regions to process
        wordRes selectRegions_;

        //- Names of fields to process
        wordRes selectFields_;

        //- Decompose polyhedra
        bool decompose_;

        //- Pointers to the requested mesh regions
        HashTable<const fvMesh*> meshes_;

        //- Backends for OpenFOAM to VTK translation (with internal caching)
        HashPtrTable<vtk::fvMeshAdaptor> backends_;

        //- Use channel enumeration from vtk::fvMeshAdaptor
        using channelEnum = vtk::fvMeshAdaptor::channel;


    // Protected Member Functions

        //- Return full channel name
        fileName channelName(channelEnum chan) const;

        //- Update/synchronize internals with catalyst backends
        void update();

        //- No copy construct
        fvMeshInput(const fvMeshInput&) = delete;

        //- No copy assignment
        void operator=(const fvMeshInput&) = delete;

public:

    //- Runtime type information
    ClassName("catalyst::fvMesh");


    // Constructors

        //- Construct from Time and dictionary
        fvMeshInput
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fvMeshInput() = default;


    // Member Functions

        //- Read the specification
        virtual bool read(const dictionary& dict);

        //- Update for changes of mesh or mesh point-motion
        virtual void update(polyMesh::readUpdateState state);

        //- Add available channels (with fields) to data query
        virtual label addChannels(dataQuery& dataq);

        //- Convert channels to vtkMultiBlockDataSet outputs
        virtual bool convert(dataQuery& dataq, outputChannels& outputs);

        //- Print information
        virtual Ostream& print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace catalyst
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
