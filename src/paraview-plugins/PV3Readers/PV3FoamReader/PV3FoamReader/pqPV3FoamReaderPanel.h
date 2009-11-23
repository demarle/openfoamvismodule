/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    pqPV3FoamReaderPanel

Description
    GUI modifications for the ParaView reader panel

    A custom panel for the PV3FoamReader.

SourceFiles
    pqPV3FoamReaderPanel.cxx

\*---------------------------------------------------------------------------*/
#ifndef pqPV3FoamReaderPanel_h
#define pqPV3FoamReaderPanel_h


#include "pqAutoGeneratedObjectPanel.h"

// Forward declaration of QT classes

class QCheckBox;
class QLineEdit;
class QTimer;
class QToolButton;

// Forward declaration of ParaView classes
class vtkSMSourceProxy;


/*---------------------------------------------------------------------------*\
                  Class pqPV3FoamReaderPanel Declaration
\*---------------------------------------------------------------------------*/

class pqPV3FoamReaderPanel
:
    public pqAutoGeneratedObjectPanel
{
    // Private data
    Q_OBJECT;
    typedef pqAutoGeneratedObjectPanel Superclass;

    //- ZeroTime checkbox
    QCheckBox* ZeroTime_;

    //- CacheMesh checkbox
    QCheckBox* CacheMesh_;

    //- Show Patch Names checkbox
    QCheckBox* ShowPatchNames_;

    //- IncludeSets checkbox
    QCheckBox* IncludeSets_;

    //- IncludeZones checkbox
    QCheckBox* IncludeZones_;

protected slots:

    void CacheMeshToggled();
    void ZeroTimeToggled();
    void RefreshPressed();
    void ShowPatchNamesToggled();
    void IncludeSetsToggled();
    void IncludeZonesToggled();


public:

    // Constructors

        //- Construct from components
        pqPV3FoamReaderPanel(pqProxy*, QWidget*);


    //- Destructor
    // virtual ~pqPV3FoamReaderPanel();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
