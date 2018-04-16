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
#include "cellCellStencilObject.H"

// VTK includes
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::fvMeshAdaptor::convertGeometryInternal()
{
    // MESH = "internal"
    if (!usingVolume())
    {
        cachedVtu_.clear();
        return;
    }

    const auto& longName = internalName;

    foamVtuData& vtuData = cachedVtu_(longName);

    vtkSmartPointer<vtkUnstructuredGrid> vtkgeom;
    if (vtuData.nPoints())
    {
        if (meshState_ == polyMesh::UNCHANGED)
        {
            if (debug)
            {
                Info<< "reuse " << longName << nl;
            }
            vtuData.reuse(); // reuse
            return;
        }
        else if (meshState_ == polyMesh::POINTS_MOVED)
        {
            if (debug)
            {
                Info<< "move points " << longName << nl;
            }
            vtkgeom = vtuData.getCopy();
            vtkgeom->SetPoints(vtuData.points(mesh_));
        }
    }

    if (!vtkgeom)
    {
        if (debug)
        {
            Info<< "Nothing usable from cache - create new geometry" << nl;
        }

        // Nothing usable from cache - create new geometry
        vtkgeom = vtuData.internal(mesh_, decomposePoly_);
    }

    vtuData.set(vtkgeom);
}


void Foam::vtk::fvMeshAdaptor::convertGeometryPatches()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const label npatches = this->nPatches();

    for (label patchId=0; patchId < npatches; ++patchId)
    {
        const polyPatch& pp = patches[patchId];
        const word& longName = pp.name();

        foamVtpData& vtpData = cachedVtp_(longName);

        vtkSmartPointer<vtkPolyData> vtkgeom;
        if (vtpData.nPoints())
        {
            if (meshState_ == polyMesh::UNCHANGED)
            {
                // Without movement is easy.
                if (debug)
                {
                    Info<< "reuse " << longName << nl;
                }
                vtpData.reuse();
                continue;
            }
            else if (meshState_ == polyMesh::POINTS_MOVED)
            {
                // Point movement on single patch is OK

                const labelList& patchIds = vtpData.additionalIds();
                if (patchIds.size() == 1)
                {
                    vtkgeom = vtpData.getCopy();
                    vtkgeom->SetPoints
                    (
                        vtk::Tools::Patch::points(patches[patchIds[0]])
                    );
                    continue;
                }
            }
        }

        vtpData.clear(); // Remove any old mappings

        if (debug)
        {
            Info<< "Creating VTK mesh for patch [" << patchId <<"] "
                << longName << endl;
        }

        // Store good patch id as additionalIds
        vtpData.additionalIds() = {patchId};

        // This is somewhat inconsistent, since we currently only have
        // normal (non-grouped) patches but this may change in the future.

        vtkgeom = vtk::Tools::Patch::mesh(patches[patchId]);

        if (vtkgeom)
        {
            vtpData.set(vtkgeom);
        }
        else
        {
            // Catch any problems
            cachedVtp_.erase(longName);
        }
    }
}


// These need to be rebuild here, since the component mesh may just have point
// motion without topology changes.
void Foam::vtk::fvMeshAdaptor::applyGhosting()
{
    if (!usingVolume())
    {
        return;
    }

    const auto* stencilPtr =
        mesh_.lookupObjectPtr<cellCellStencilObject>
        (
            cellCellStencilObject::typeName
        );

    if (!stencilPtr)
    {
        return;
    }

    const auto& longName = internalName;

    auto iter = cachedVtu_.find(longName);
    if (!iter.found() || !iter.object().dataset)
    {
        // Should not happen, but for safety require a vtk geometry
        return;
    }
    foamVtuData& vtuData = iter.object();
    auto dataset = vtuData.dataset;

    const auto& stencil = *stencilPtr;
    const labelUList& cellMap = vtuData.cellMap();

    auto vtkgcell = dataset->GetCellGhostArray();
    if (!vtkgcell)
    {
        vtkgcell = dataset->AllocateCellGhostArray();
    }

    // auto vtkgpoint = dataset->GetPointGhostArray();
    // if (!vtkgpoint)
    // {
    //     vtkgpoint = dataset->AllocatePointGhostArray();
    // }

    UList<uint8_t> gcell =
        vtk::Tools::asUList(vtkgcell, cellMap.size());

    vtkgcell->FillValue(0); // Initialize to zero

    const labelUList& types = stencil.cellTypes();

    forAll(cellMap, i)
    {
        const label cellType = types[cellMap[i]];

        if (cellType == cellCellStencil::INTERPOLATED)
        {
            gcell[i] |=
            (
                vtkDataSetAttributes::DUPLICATECELL
            );
        }
        else if (cellType == cellCellStencil::HOLE)
        {
            // Need duplicate (not just HIDDENCELL) for it to be properly
            // recognized
            gcell[i] |=
            (
                vtkDataSetAttributes::DUPLICATECELL
              | vtkDataSetAttributes::HIDDENCELL
            );
        }
    }

    dataset->GetCellData()->AddArray(vtkgcell);
}


// ************************************************************************* //
