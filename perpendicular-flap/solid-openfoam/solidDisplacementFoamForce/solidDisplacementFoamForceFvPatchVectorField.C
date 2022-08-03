#include "solidDisplacementFoamForceFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidDisplacementFoamForceFvPatchVectorField::
solidDisplacementFoamForceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    force_(p.size(), vector::zero),
    forceFieldPtr_(),
    curTimeIndex_(-1)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


solidDisplacementFoamForceFvPatchVectorField::
solidDisplacementFoamForceFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    force_(p.size(), vector::zero),
    forceFieldPtr_(),
    curTimeIndex_(-1)
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    // Initialise traction and pressure to zero
    traction() = vector::zero;
    pressure() = 0.0;

    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    // Check how force is defined
    if (dict.found("force") && dict.found("forceField"))
    {
        FatalErrorIn
        (
            "solidDisplacementFoamForceFvPatchVectorField::solidDisplacementFoamForceFvPatchVectorField"
        )   << "Only force or forceField can be "
            << "specified, not both!"
            << abort(FatalError);
    }
    else if (dict.found("forceField"))
    {
        Info<< "    force is specified as a field" << endl;
        forceFieldPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    word(dict.lookup("forceField")),
                    patch().boundaryMesh().mesh().time().timeName(),
                    patch().boundaryMesh().mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                patch().boundaryMesh().mesh()
            )
        );
    }
    else
    {
        force_ = vectorField("force", dict, p.size());
    }
}


solidDisplacementFoamForceFvPatchVectorField::
solidDisplacementFoamForceFvPatchVectorField
(
    const solidDisplacementFoamForceFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(stpvf, p, iF, mapper),
#ifdef OPENFOAMFOUNDATION
    force_(mapper(stpvf.force_)),
#else
    force_(stpvf.force_, mapper),
#endif
    forceFieldPtr_(),
    curTimeIndex_(stpvf.curTimeIndex_)
{}

#ifndef OPENFOAMFOUNDATION
solidDisplacementFoamForceFvPatchVectorField::solidDisplacementFoamForceFvPatchVectorField
(
    const solidDisplacementFoamForceFvPatchVectorField& stpvf
)
:
    tractionDisplacementFvPatchVectorField(stpvf),
    force_(stpvf.force_),
    forceFieldPtr_(),
    curTimeIndex_(stpvf.curTimeIndex_)
{}
#endif

solidDisplacementFoamForceFvPatchVectorField::solidDisplacementFoamForceFvPatchVectorField
(
    const solidDisplacementFoamForceFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(stpvf, iF),
    force_(stpvf.force_),
    forceFieldPtr_(),
    curTimeIndex_(stpvf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidDisplacementFoamForceFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);

#ifdef OPENFOAMFOUNDATION
    m(force_, force_);
#else
    force_.autoMap(m);
#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidDisplacementFoamForceFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);

    const solidDisplacementFoamForceFvPatchVectorField& dmptf =
        refCast<const solidDisplacementFoamForceFvPatchVectorField>(ptf);

    force_.rmap(dmptf.force_, addr);
}


// Update the coefficients associated with the patch field
void solidDisplacementFoamForceFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != db().time().timeIndex())
    {
        curTimeIndex_ = db().time().timeIndex();

        // Called once per time-step
        if (forceFieldPtr_.valid())
        {
            // Force the traction field boundary conditions to update
            const_cast<volVectorField&>
            (
                forceFieldPtr_()
            ).correctBoundaryConditions();
        }
    }

    if (forceFieldPtr_.valid())
    {
        force_ = forceFieldPtr_().boundaryField()[patch().index()];
    }
    
    // Convert the force field to a traction field
    // Note: this assumes small strains / linear geometry
    traction() = force_/patch().magSf();

    // Apply traction
    tractionDisplacementFvPatchVectorField::updateCoeffs();
}

void solidDisplacementFoamForceFvPatchVectorField::write(Ostream& os) const
{
    // Bug-fix: courtesy of Michael@UW at https://www.cfd-online.com/Forums/
    // openfoam-cc-toolkits-fluid-structure-interaction/221892-solved-paraview
    // -cant-read-solids-files-duplicate-entries-keyword-value.html#post762325
    //tractionDisplacementFvPatchVectorField::write(os);
    fvPatchVectorField::write(os);

    if (forceFieldPtr_.valid())
    {
        os.writeKeyword("forceField")
            << forceFieldPtr_().name() << token::END_STATEMENT << nl;
    }
    else
    {
#ifdef OPENFOAMFOUNDATION
        writeEntry(os, "force", force_);
#else
        force_.writeEntry("force", os);
#endif
    }

#ifdef OPENFOAMFOUNDATION
    writeEntry(os, "value", *this);
    writeEntry(os, "gradient", gradient());
#else
    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidDisplacementFoamForceFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
