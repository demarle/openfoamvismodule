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
\*---------------------------------------------------------------------------*/

#include "catalystCoprocess.H"
#include "Time.H"
#include "stringOps.H"
#include "OSspecific.H"

// VTK includes
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkMultiBlockDataSet.h>

// Catalyst includes
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(catalystCoprocess, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
namespace
{

#if 0
static void printInfo(vtkCPDataDescription* descrip)
{
    if (!descrip)
    {
        return;
    }

    const unsigned nItems = descrip->GetNumberOfInputDescriptions();

    for (unsigned itemi = 0; itemi < nItems; ++itemi)
    {
        vtkCPInputDataDescription* input = descrip->GetInputDescription(itemi);
        if (!input) continue;  // should not happen

        Info<<"input: " << descrip->GetInputDescriptionName(itemi) << nl;

        const unsigned nFields = input->GetNumberOfFields();
        for (unsigned fieldi = 0; fieldi < nFields; ++fieldi)
        {
            Info<< "    field: " << input->GetFieldName(fieldi) << nl;
        }
        if (!nFields) Info<<"     no fields requested" << nl;
    }
}
#endif

} // End anonymous namespace
} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::catalystCoprocess::expand
(
    List<string>& scripts,
    const dictionary& dict
)
{
    label nscript = 0;

    forAll(scripts, scripti)
    {
        string& s = scripts[scripti];

        stringOps::inplaceExpand(s, dict, true, true);
        fileName::clean(s);  // Remove trailing, repeated slashes etc.

        if (isFile(s))
        {
            if (nscript != scripti)
            {
                scripts[nscript] = std::move(s);
            }
            ++nscript;
        }
    }

    scripts.resize(nscript);

    return nscript;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class DataType>
bool Foam::catalystCoprocess::processImpl
(
    const dataQuery& dataq,
    vtkSmartPointer<DataType>& output
)
{
    vtkCPDataDescription* descrip = dataq.get();

    if (!coproc_->RequestDataDescription(descrip))
    {
        return false;
    }

    for (const word& chanName : dataq.channels())
    {
        auto* input = descrip->GetInputDescriptionByName(chanName.c_str());

        if (input && input->GetIfGridIsNecessary())
        {
            input->SetGrid(output);
        }
    }

    coproc_->CoProcess(descrip);
    return true;
}


template<class DataType>
bool Foam::catalystCoprocess::processImpl
(
    const dataQuery& dataq,
    HashTable<vtkSmartPointer<DataType>>& outputs
)
{
    vtkCPDataDescription* descrip = dataq.get();

    if (!coproc_->RequestDataDescription(descrip))
    {
        return false;
    }

    for (const word& chanName : dataq.channels())
    {
        if (outputs.found(chanName))
        {
            auto* input = descrip->GetInputDescriptionByName(chanName.c_str());

            if (input && input->GetIfGridIsNecessary())
            {
                input->SetGrid(outputs[chanName]);
            }
        }
    }

    coproc_->CoProcess(descrip);
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::catalystCoprocess::timeQuery::timeQuery(const Foam::Time& t)
:
    timeValue(t.timeOutputValue()),
    timeIndex(t.timeIndex()),
    forced(t.timeOutputValue() >= t.endTime().value())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::catalystCoprocess::~catalystCoprocess()
{
    // stop(), but without output
    if (coproc_)
    {
        coproc_->Delete();
        coproc_ = nullptr;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::catalystCoprocess::good() const
{
    return coproc_;
}


void Foam::catalystCoprocess::stop()
{
    if (coproc_)
    {
        coproc_->Delete();
        coproc_ = nullptr;
        Info<< "Stop catalyst" << endl;
    }
}


void Foam::catalystCoprocess::reset(const fileName& outputDir)
{
    #ifdef USE_CATALYST_WORKING_DIRECTORY
    if (coproc_ == nullptr)
    {
        coproc_ = vtkCPProcessor::New();
        coproc_->Initialize(outputDir.c_str());
        Info<< "Connecting ParaView Catalyst..." << endl;
    }
    else
    {
        coproc_->RemoveAllPipelines();

        if (outputDir == coproc_->GetWorkingDirectory())
        {
            Info<< "Rebinding ParaView Catalyst..." << endl;
        }
        else
        {
            // Changed working directory ... redo everything.
            coproc_->Delete();
            coproc_ = nullptr;

            reset(outputDir);
        }
    }
    #else
    if (coproc_ == nullptr)
    {
        coproc_ = vtkCPProcessor::New();
        coproc_->Initialize();
    }
    else
    {
        coproc_->RemoveAllPipelines();
        Info<< "Rebinding ParaView Catalyst..." << endl;
    }
    Info<< "    Caution: using current working directory" << nl
        << "    which may not be the same as the simulation directory" << endl;
    #endif
}


void Foam::catalystCoprocess::reset
(
    const fileName& outputDir,
    const UList<string>& scripts
)
{
    reset(outputDir);

    int nscript = 0;
    for (const auto& script : scripts)
    {
        Info<< "Adding pipeline[" << nscript << "] " << script << endl;
        ++nscript;

        vtkNew<vtkCPPythonScriptPipeline> pipeline;
        pipeline->Initialize(script.c_str());
        coproc_->AddPipeline(pipeline.GetPointer());
    }

    // Do something different with (nscript == 0) ?
}


Foam::HashTable<Foam::wordHashSet>
Foam::catalystCoprocess::query
(
    dataQuery& dataq,
    const wordHashSet& allFields
)
{
    // Desirable to also know which fields have been requested
    HashTable<wordHashSet> requests;

    if (!good())
    {
        Info<< "No ParaView Catalyst initialized" << endl;
        return requests;
    }

    if (dataq.channels().empty())
    {
        // No channels names have been published by the simulation
        return requests;
    }

    vtkCPDataDescription* descrip = dataq.get();

    descrip->SetTimeData(dataq.timeValue, dataq.timeIndex);
    descrip->SetForceOutput(dataq.forced);

    // Sort out which channels already exist, are new, or disappeared
    {
        // The currently defined channels
        wordHashSet currChannels;

        const unsigned n = descrip->GetNumberOfInputDescriptions();
        for (unsigned i=0; i < n; ++i)
        {
            currChannels.insert
            (
                word::validate(descrip->GetInputDescriptionName(i))
            );
        }

        wordHashSet newChannels(dataq.channels());
        wordHashSet oldChannels(currChannels);
        oldChannels.erase(newChannels);

        if (oldChannels.size())
        {
            descrip->ResetAll();
        }
        else
        {
            newChannels.erase(currChannels);
        }

        // Add channels
        for (const word& chanName : newChannels)
        {
            descrip->AddInput(chanName.c_str());
            auto* input = descrip->GetInputDescriptionByName(chanName.c_str());

            for (const word& fieldName : allFields)
            {
                input->AddPointField(fieldName.c_str());
                input->AddCellField(fieldName.c_str());
            }
        }

        // Note: this misses updating field information for previously
        // existing inputs.
    }

    if
    (
        !coproc_->RequestDataDescription(descrip)
     || !descrip->GetIfAnyGridNecessary()
    )
    {
        return requests;
    }

    for (const word& chanName : dataq.channels())
    {
        auto* input = descrip->GetInputDescriptionByName(chanName.c_str());

        if (input && input->GetIfGridIsNecessary())
        {
            wordHashSet& fields = requests(chanName);  // auto-vivify

            for (const word& fieldName : allFields)
            {
                if (input->IsFieldNeeded(fieldName.c_str()))
                {
                    fields.insert(fieldName);
                }
            }
        }
    }

    return requests;
}


bool Foam::catalystCoprocess::process
(
    const dataQuery& dataq,
    vtkSmartPointer<vtkMultiBlockDataSet>& output
)
{
    return processImpl(dataq, output);
}


bool Foam::catalystCoprocess::process
(
    const dataQuery& dataq,
    HashTable<vtkSmartPointer<vtkMultiBlockDataSet>>& outputs
)
{
    return processImpl(dataq, outputs);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const catalystCoprocess::timeQuery& when
)
{
    os << "Time = " << when.timeValue << ", index: " << when.timeIndex;
    return os;
}


// ************************************************************************* //
