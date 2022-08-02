/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    solidDisplacementFoamForceFvPatchVectorField

Description
    Force boundary condition for solidDisplacementFoam solvers, ported from
    the solidForce condition in solids4foam.

    This is similar to tractionDisplacement except the absolute force per face
    is specified.

    Note: this applies the force per face, not to the entire patch.

    Obviously, for a uniform force field, the total force applied to the patch
    as the mesh is refined as the force per face stays constant.

    The primary motivation for this condition is for use FSI cases using the
    preCICE coupling library.

    The force field can be directly specified or specified via a "force"
    volVectorField, as used with preCICE.

Author
    Philip Cardiff, UCD.

SourceFiles
    solidDisplacementFoamForceFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidDisplacementFoamForceFvPatchVectorField_H
#define solidDisplacementFoamForceFvPatchVectorField_H

#include "tractionDisplacementFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class solidDisplacementFoamForceFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidDisplacementFoamForceFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{

    // Private Data

        //- Force
        vectorField force_;

        //- Force field given as a vol field
        autoPtr<volVectorField> forceFieldPtr_;

        //- Current time index
        label curTimeIndex_;

public:

    //- Runtime type information
    TypeName("solidDisplacementFoamForce");


    // Constructors

        //- Construct from patch and internal field
        solidDisplacementFoamForceFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidDisplacementFoamForceFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidDisplacementFoamForceFvPatchVectorField onto a new patch
        solidDisplacementFoamForceFvPatchVectorField
        (
            const solidDisplacementFoamForceFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidDisplacementFoamForceFvPatchVectorField
        (
            const solidDisplacementFoamForceFvPatchVectorField&
#ifdef OPENFOAMFOUNDATION
        ) = delete;
#else
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new solidDisplacementFoamForceFvPatchVectorField(*this)
            );
        }
#endif

        //- Construct as copy setting internal field reference
        solidDisplacementFoamForceFvPatchVectorField
        (
            const solidDisplacementFoamForceFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new solidDisplacementFoamForceFvPatchVectorField(*this, iF)
            );
        }



    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
