/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "foamVtkFvMeshAdaptor.H"
#include "Pstream.H"

// OpenFOAM includes
#include "fvMesh.H"

// VTK includes
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkInformation.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{
    defineTypeNameAndDebug(fvMeshAdaptor, 0);
}
}

const Foam::Enum
<
    Foam::vtk::fvMeshAdaptor::channel
>
Foam::vtk::fvMeshAdaptor::channelNames
{
    { channel::MESH,    "mesh" },
    { channel::PATCHES, "patches" },
    { channel::INPUT,   "input" },
};


const Foam::word Foam::vtk::fvMeshAdaptor::internalName("internal");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::fvMeshAdaptor::fvMeshAdaptor
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    channels_(INPUT),
    decomposePoly_(false),
    meshState_(polyMesh::TOPO_CHANGE)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::fvMeshAdaptor::channels(const wordList& chanNames)
{
    unsigned chanIds = 0;
    for (const word& chan : chanNames)
    {
        if (channelNames.found(chan))
        {
            chanIds |= channelNames[chan];
        }
    }

    channels(chanIds);
}


void Foam::vtk::fvMeshAdaptor::channels(unsigned chanIds)
{
    channels_ = (chanIds & 0x3);

    if (!usingVolume())
    {
        cachedVtu_.clear();
    }

    if (!usingPatches())
    {
        cachedVtp_.clear();
    }
}


Foam::label Foam::vtk::fvMeshAdaptor::channels() const
{
    return label(channels_);
}


bool Foam::vtk::fvMeshAdaptor::usingVolume() const
{
    // MESH = "internal"
    return (channels_ & (INPUT | MESH));
}


bool Foam::vtk::fvMeshAdaptor::usingPatches() const
{
    // PATCHES
    return (channels_ & (INPUT | PATCHES));
}



Foam::label Foam::vtk::fvMeshAdaptor::nPatches() const
{
    // Restrict to non-processor patches.
    // This value is invariant across all processors.

    if (usingPatches())
    {
        return mesh_.boundaryMesh().nNonProcessor();
    }

    return 0;
}


void Foam::vtk::fvMeshAdaptor::updateContent(const wordRes& selectFields)
{
    const bool oldDecomp = decomposePoly_;

    // TODO from dictionary
    //decomposePoly_ = !reader_->GetUseVTKPolyhedron();

    // Update cached, saved, unneed values.

    HashSet<string> nowActive;

    // MESH = "internal"
    if (usingVolume())
    {
        nowActive.insert(internalName);
    }

    // PATCHES
    // Restrict to non-processor patches.
    // This value is invariant across all processors.

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const label npatches = this->nPatches();

    for (label patchId=0; patchId < npatches; ++patchId)
    {
        const polyPatch& pp = patches[patchId];
        nowActive.insert(pp.name());
    }

    // Dispose of unneeded components
    cachedVtp_.retain(nowActive);
    cachedVtu_.retain(nowActive);

    if
    (
        meshState_ == polyMesh::TOPO_CHANGE
     || meshState_ == polyMesh::TOPO_PATCH_CHANGE
    )
    {
        // Eliminate cached values that would be unreliable
        forAllIters(cachedVtp_, iter)
        {
            iter.object().clearGeom();
            iter.object().clear();
        }
        forAllIters(cachedVtu_, iter)
        {
            iter.object().clearGeom();
            iter.object().clear();
        }
    }
    else if (oldDecomp != decomposePoly_)
    {
        // poly-decompose changed - dispose of cached values
        forAllIters(cachedVtu_, iter)
        {
            iter.object().clearGeom();
            iter.object().clear();
        }
    }

    convertGeometryInternal();
    convertGeometryPatches();
    applyGhosting();
    convertVolFields(selectFields);
    meshState_ = polyMesh::UNCHANGED;
}


vtkSmartPointer<vtkMultiBlockDataSet>
Foam::vtk::fvMeshAdaptor::output(const wordRes& select)
{
    updateContent(select);

    // All individual datasets are vtkMultiPieceDataSet for improved
    // handling downstream.

    label rank = 0;
    label nproc = 1;

    if (Pstream::parRun())
    {
        rank  = Pstream::myProcNo();
        nproc = Pstream::nProcs();
    }


    auto outputs = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    // MESH
    {
        unsigned int blockNo = 0;

        do
        {
            const auto& longName = internalName;
            auto iter = cachedVtu_.find(longName);
            if (!iter.found() || !iter.object().dataset)
            {
                Pout<<"Cache miss for VTU " << longName << endl;
                break; // Should never happen
            }

            foamVtuData& vtuData = iter.object();

            auto pieces = vtkSmartPointer<vtkMultiPieceDataSet>::New();

            pieces->SetNumberOfPieces(nproc);
            pieces->SetPiece(rank, vtuData.dataset);

            outputs->SetBlock(blockNo, pieces);

            outputs->GetMetaData(blockNo)->Set
            (
                vtkCompositeDataSet::NAME(),
                internalName
            );

            ++blockNo;
        } while (false);  // do once
    }

    // PATCHES
    const label npatches = this->nPatches();
    if (npatches)
    {
        unsigned int blockNo = 0;

        auto output = vtkSmartPointer<vtkMultiBlockDataSet>::New();

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        for (label patchId=0; patchId < npatches; ++patchId)
        {
            const polyPatch& pp = patches[patchId];
            const word& longName = pp.name();

            auto iter = cachedVtp_.find(longName);
            if (!iter.found() || !iter.object().dataset)
            {
                Pout<<"Cache miss for VTP patch " << longName << endl;
                break; // Should never happen
            }

            foamVtpData& vtpData = iter.object();

            auto pieces = vtkSmartPointer<vtkMultiPieceDataSet>::New();

            pieces->SetNumberOfPieces(nproc);
            pieces->SetPiece(rank, vtpData.dataset);

            output->SetBlock(blockNo, pieces);

            output->GetMetaData(blockNo)->Set
            (
                vtkCompositeDataSet::NAME(),
                longName
            );

            ++blockNo;
        }

        outputs->SetBlock(1, output);
        outputs->GetMetaData(1)->Set
        (
            vtkCompositeDataSet::NAME(),
            channelNames[channel::PATCHES]
        );
    }

    // Would actually like to have this:
    //     outputs->SetName(mesh_.name().c_str());
    // but do that in the caller side

    return outputs;
}


void Foam::vtk::fvMeshAdaptor::updateState(polyMesh::readUpdateState state)
{
    // Only move to worse states
    switch (state)
    {
        case polyMesh::UNCHANGED:
            break;

        case polyMesh::POINTS_MOVED:
            if (meshState_ == polyMesh::UNCHANGED)
            {
                meshState_ = polyMesh::POINTS_MOVED;
            }
            break;

        case polyMesh::TOPO_CHANGE:
        case polyMesh::TOPO_PATCH_CHANGE:
            meshState_ = polyMesh::TOPO_CHANGE;
            break;
    }
}


// ************************************************************************* //
