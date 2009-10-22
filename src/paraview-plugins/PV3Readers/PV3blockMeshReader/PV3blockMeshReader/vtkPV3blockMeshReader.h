/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    vtkPV3blockMeshReader

Description
    reads a dataset in OpenFOAM bockMesh format

    vtkPV3blockMeshReader creates an multiblock dataset.
    It uses the OpenFOAM infrastructure (blockMesh).

SourceFiles
    vtkPV3blockMeshReader.cxx

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3blockMeshReader_h
#define vtkPV3blockMeshReader_h

// VTK includes
#include "vtkMultiBlockDataSetAlgorithm.h"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

// VTK forward declarations
class vtkDataArraySelection;
class vtkCallbackCommand;

namespace Foam
{
    class vtkPV3blockMesh;
}

/*---------------------------------------------------------------------------*\
                   Class vtkPV3blockMeshReader Declaration
\*---------------------------------------------------------------------------*/

class VTK_IO_EXPORT vtkPV3blockMeshReader
:
    public vtkMultiBlockDataSetAlgorithm
{
public:
    vtkTypeRevisionMacro(vtkPV3blockMeshReader,vtkMultiBlockDataSetAlgorithm);
    void PrintSelf(ostream&, vtkIndent);

    static vtkPV3blockMeshReader* New();

    // Description:
    // Set/Get the filename.
    vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    // Description:
    // Display corner point labels
    virtual void SetShowPointNumbers(int);
    vtkGetMacro(ShowPointNumbers, int);

    // Description:
    // GUI update control
    vtkSetMacro(UpdateGUI, int);
    vtkGetMacro(UpdateGUI, int);


    // Description:
    // Parts (blocks) selection list control
    vtkDataArraySelection* GetPartSelection();
    int  GetNumberOfPartArrays();
    int  GetPartArrayStatus(const char*);
    void SetPartArrayStatus(const char*, int status);
    const char* GetPartArrayName(int index);

    // Description:
    // Parts (blocks) selection list control
    vtkDataArraySelection* GetCurvedEdgesSelection();
    int  GetNumberOfCurvedEdgesArrays();
    int  GetCurvedEdgesArrayStatus(const char*);
    void SetCurvedEdgesArrayStatus(const char*, int status);
    const char* GetCurvedEdgesArrayName(int index);

    // Description:
    // Callback registered with the SelectionObserver
    // for all the selection lists
    static void SelectionModifiedCallback
    (
        vtkObject* caller,
        unsigned long eid,
        void* clientdata,
        void* calldata
    );

    void SelectionModified();


protected:

    //- Construct null
    vtkPV3blockMeshReader();

    //- Destructor
    ~vtkPV3blockMeshReader();

    //- Return information about mesh, times, etc without loading anything
    virtual int RequestInformation
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    );

    //- Get the mesh/fields for a particular time
    virtual int RequestData
    (
        vtkInformation*,
        vtkInformationVector**,
        vtkInformationVector*
    );

    //- Fill in additional port information
    virtual int FillOutputPortInformation(int, vtkInformation*);

    // The observer to modify this object when array selections are modified
    vtkCallbackCommand* SelectionObserver;

    char* FileName;


private:

    //- Disallow default bitwise copy construct
    vtkPV3blockMeshReader(const vtkPV3blockMeshReader&);

    //- Disallow default bitwise assignment
    void operator=(const vtkPV3blockMeshReader&);

    //- Add/remove point numbers to/from the view
    void updatePointNumbersView(const bool show);


    //- Show Point Numbers
    int ShowPointNumbers;

    //- Dummy variable/switch to invoke a reader update
    int UpdateGUI;

    vtkDataArraySelection* PartSelection;

    vtkDataArraySelection* CurvedEdgesSelection;

    //BTX
    Foam::vtkPV3blockMesh* foamData_;
    //ETX
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
