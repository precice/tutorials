/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

#include "codedFvOptionTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = fa508e6ae0244268d3eb64a1b0f1e3afb6179e04
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void sourceTime_fa508e6ae0244268d3eb64a1b0f1e3afb6179e04(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(sourceTimeFvOptionvectorSource, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    sourceTimeFvOptionvectorSource,
    dictionary
);

} // End namespace fv
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::
sourceTimeFvOptionvectorSource::
sourceTimeFvOptionvectorSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh)
{
    if (false)
    {
        printMessage("Construct sourceTime fvOption from dictionary");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::
sourceTimeFvOptionvectorSource::
~sourceTimeFvOptionvectorSource()
{
    if (false)
    {
        printMessage("Destroy sourceTime");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::fv::
sourceTimeFvOptionvectorSource::correct
(
    GeometricField<vector, fvPatchField, volMesh>& fld
)
{
    if (false)
    {
        Info<< "sourceTimeFvOptionvectorSource::correct()\n";
    }

//{{{ begin code
    #line 25 "/home/tirgendetwas/preCICE/tutorials/volume-coupled-flow/fluid-openfoam/constant/fvOptions.codedSource"
const labelList& cells = this->cells();
        const volVectorField& U_vol = mesh_.lookupObject<volVectorField>("U_vol");
        for(auto cell : cells)
        {
            fld[cell].x() = U_vol[cell].x();
        }
//}}} end code
}


void
Foam::fv::
sourceTimeFvOptionvectorSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<< "sourceTimeFvOptionvectorSource::addSup()\n";
    }

//{{{ begin code - warn/fatal if not implemented?
    #line 35 "/home/tirgendetwas/preCICE/tutorials/volume-coupled-flow/fluid-openfoam/constant/fvOptions.codedSource"
return;
//}}} end code
}


void
Foam::fv::
sourceTimeFvOptionvectorSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<< "sourceTimeFvOptionvectorSource::addSup(rho)\n";
    }

//{{{ begin code - warn/fatal if not implemented?
    #line 40 "/home/tirgendetwas/preCICE/tutorials/volume-coupled-flow/fluid-openfoam/constant/fvOptions.codedSource"
return;
//}}} end code
}


void
Foam::fv::
sourceTimeFvOptionvectorSource::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<< "sourceTimeFvOptionvectorSource::constrain()\n";
    }

//{{{ begin code
    #line 20 "/home/tirgendetwas/preCICE/tutorials/volume-coupled-flow/fluid-openfoam/constant/fvOptions.codedSource"
return;
//}}} end code
}


// ************************************************************************* //

