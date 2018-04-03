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
    Foam::catalystCoprocess

Description
    Low-level interface between OpenFOAM and ParaView Catalyst.

    \code
    ... initialize catalyst

    ... define a data query for catalyst

    // Data description for co-processing
    vtkNew<vtkCPDataDescription> descrip;

    // Form data query for catalyst
    catalystCoprocess::dataQuery dataq(channelNames, runTime, descrip.Get());

    // Query catalyst
    HashTable<wordHashSet> expecting = adaptor_().query(dataq, fields);

    \endcode

SourceFiles
    catalystCoprocess.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_catalystCoprocess_H
#define functionObjects_catalystCoprocess_H

#include "className.H"
#include "wordList.H"
#include "stringList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward Declarations
class vtkCPProcessor;
class vtkCPDataDescription;
class vtkMultiBlockDataSet;
template<class T> class vtkSmartPointer;

namespace Foam
{

// Forward Declarations
class dictionary;
class Time;
class Ostream;
class catalystCoprocess;

/*---------------------------------------------------------------------------*\
                      Class catalystCoprocess Declaration
\*---------------------------------------------------------------------------*/

class catalystCoprocess
{
public:

    //- Simple structure for time queries
    struct timeQuery
    {
        double timeValue;
        long   timeIndex;
        bool   forced;

        constexpr timeQuery
        (
            double val,
            long index,
            bool forcedOutput=false
        ) noexcept
        :
            timeValue(val),
            timeIndex(index),
            forced(forcedOutput)
        {}

        timeQuery(const Foam::Time& currTime);
    };


    //- Simple structure for data description queries
    //  The storage for the description is held outside this class
    class dataQuery
    :
        public timeQuery
    {
        //- Catalyst channel names to query
        List<word> channels_;
        mutable vtkCPDataDescription* descrip_;

    public:

        dataQuery
        (
            const UList<word>& channelNames,
            const timeQuery& when,
            vtkCPDataDescription* description
        )
        :
            timeQuery(when),
            channels_(channelNames),
            descrip_(description)
        {}

        dataQuery
        (
            const UList<word>& channelNames,
            const Foam::Time& when,
            vtkCPDataDescription* description
        )
        :
            dataQuery(channelNames, timeQuery(when), description)
        {}

        //- The declared output channels
        const UList<word>& channels() const
        {
            return channels_;
        }

        //- Pointer to coprocess data description
        vtkCPDataDescription* get() const
        {
            return this->descrip_;
        }
    };


private:

    // Private Data

        //- The currently active ParaView Catalyst co-process
        vtkCPProcessor* coproc_;


    // Private Member Functions

        //- Process single output channel
        template<class DataType>
        bool processImpl
        (
            const dataQuery& query,
            vtkSmartPointer<DataType>& outputs
        );

        //- Process multiple output channels
        template<class DataType>
        bool processImpl
        (
            const dataQuery& query,
            HashTable<vtkSmartPointer<DataType>>& outputs
        );

public:

    // Static Methods

        //- Expand strings as filenames, retaining only those that exist
        static label expand(List<string>& scripts, const dictionary& dict);


    //- Define class name and debug
    ClassName("catalyst");


    // Constructors

        //- Construct null. Does not initialize catalyst.
        constexpr catalystCoprocess() noexcept
        :
            coproc_(nullptr)
        {}


    //- Destructor. Shutdown process.
    ~catalystCoprocess();


    // Member Functions

        //- \return True if the coprocess has been initialized.
        bool good() const;

        //- Reset/initialize pipeline without pipeline scripts.
        void reset();

        //- Reset/initialize pipeline with python scripts.
        //  The script names must have already been resolved prior calling this.
        void reset(const UList<string>& scripts);

        //- Query the coprocess pipelines if they should be executed at this
        //- iteration and possibly which fields they require.
        //
        // \param[in,out] dataq the data query for catalyst.
        //     On input it contains the published channel names, the current
        //     simulation time (index, value) and allocation for the coprocess
        //     data description.
        //     On output the data description will be filled with the field
        //     names added per channel.
        // \param[in] allFields the fields that can be published from the
        //     simulation.
        //
        // \return HashTable with fields requested (what Catalyst expects)
        //     on a per-channel basis.
        HashTable<wordHashSet> query
        (
            dataQuery& dataq,
            const wordHashSet& allFields
        );

        //- Single-channel source (eg, "input" or "cloud", ...)
        bool process
        (
            const dataQuery& dataq,
            vtkSmartPointer<vtkMultiBlockDataSet>& output
        );

        //- Multi-channel source (eg, "input", "mesh", "patches")
        bool process
        (
            const dataQuery& dataq,
            HashTable<vtkSmartPointer<vtkMultiBlockDataSet>>& outputs
        );

        //- Finalize
        void stop();

};


// Ostream Operator
Ostream& operator<<(Ostream& os, const catalystCoprocess::timeQuery& when);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
