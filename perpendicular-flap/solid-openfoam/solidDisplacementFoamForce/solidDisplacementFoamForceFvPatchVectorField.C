#include "addToRunTimeSelectionTable.H"
#include "solidDisplacementFoamForceFvPatchVectorField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam {

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidDisplacementFoamForceFvPatchVectorField::
    solidDisplacementFoamForceFvPatchVectorField(
        const fvPatch                           &p,
        const DimensionedField<vector, volMesh> &iF)
    : fixedGradientFvPatchVectorField(p, iF),
      force_(p.size(), vector::zero),
      forceFieldPtr_(),
      curTimeIndex_(-1)
{
  fvPatchVectorField::operator=(patchInternalField());
  gradient() = vector::zero;
}

solidDisplacementFoamForceFvPatchVectorField::
    solidDisplacementFoamForceFvPatchVectorField(
        const fvPatch                           &p,
        const DimensionedField<vector, volMesh> &iF,
        const dictionary                        &dict)
    : fixedGradientFvPatchVectorField(p, iF),
      force_(p.size(), vector::zero),
      forceFieldPtr_(),
      curTimeIndex_(-1)
{
  Info << "Creating " << type() << " boundary condition" << endl;

  if (dict.found("gradient")) {
    gradient() = vectorField("gradient", dict, p.size());
  } else {
    gradient() = vector::zero;
  }

  if (dict.found("value")) {
    Field<vector>::operator=(vectorField("value", dict, p.size()));
  } else {
    fvPatchVectorField::operator=(patchInternalField());
  }

  // Check how force is defined
  if (dict.found("force") && dict.found("forceField")) {
    FatalErrorIn(
        "solidDisplacementFoamForceFvPatchVectorField::solidDisplacementFoamForceFvPatchVectorField")
        << "Only force or forceField can be "
        << "specified, not both!"
        << abort(FatalError);
  } else if (dict.found("forceField")) {
    Info << "    force is specified as a field" << endl;
    forceFieldPtr_.reset(
        new volVectorField(
            IOobject(
                word(dict.lookup("forceField")),
                patch().boundaryMesh().mesh().time().timeName(),
                patch().boundaryMesh().mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE),
            patch().boundaryMesh().mesh()));
  } else {
    force_ = vectorField("force", dict, p.size());
  }
}

solidDisplacementFoamForceFvPatchVectorField::
    solidDisplacementFoamForceFvPatchVectorField(
        const solidDisplacementFoamForceFvPatchVectorField &stpvf,
        const fvPatch                                      &p,
        const DimensionedField<vector, volMesh>            &iF,
        const fvPatchFieldMapper                           &mapper)
    : fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
#ifdef OPENFOAMFOUNDATION
      force_(mapper(stpvf.force_)),
#else
      force_(stpvf.force_, mapper),
#endif
      forceFieldPtr_(),
      curTimeIndex_(stpvf.curTimeIndex_)
{
}

#ifndef OPENFOAMFOUNDATION
solidDisplacementFoamForceFvPatchVectorField::solidDisplacementFoamForceFvPatchVectorField(
    const solidDisplacementFoamForceFvPatchVectorField &stpvf)
    : fixedGradientFvPatchVectorField(stpvf),
      force_(stpvf.force_),
      forceFieldPtr_(),
      curTimeIndex_(stpvf.curTimeIndex_)
{
}
#endif

solidDisplacementFoamForceFvPatchVectorField::solidDisplacementFoamForceFvPatchVectorField(
    const solidDisplacementFoamForceFvPatchVectorField &stpvf,
    const DimensionedField<vector, volMesh>            &iF)
    : fixedGradientFvPatchVectorField(stpvf, iF),
      force_(stpvf.force_),
      forceFieldPtr_(),
      curTimeIndex_(stpvf.curTimeIndex_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void solidDisplacementFoamForceFvPatchVectorField::autoMap(
    const fvPatchFieldMapper &m)
{
  fixedGradientFvPatchVectorField::autoMap(m);

#ifdef OPENFOAMFOUNDATION
  m(force_, force_);
#else
  force_.autoMap(m);
#endif
}

// Reverse-map the given fvPatchField onto this fvPatchField
void solidDisplacementFoamForceFvPatchVectorField::rmap(
    const fvPatchVectorField &ptf,
    const labelList          &addr)
{
  fixedGradientFvPatchVectorField::rmap(ptf, addr);

  const solidDisplacementFoamForceFvPatchVectorField &dmptf =
      refCast<const solidDisplacementFoamForceFvPatchVectorField>(ptf);

  force_.rmap(dmptf.force_, addr);
}

// Update the coefficients associated with the patch field
void solidDisplacementFoamForceFvPatchVectorField::updateCoeffs()
{
  if (updated()) {
    return;
  }

  if (curTimeIndex_ != db().time().timeIndex()) {
    curTimeIndex_ = db().time().timeIndex();

    // Called once per time-step
    if (forceFieldPtr_.valid()) {
      // Force the traction field boundary conditions to update
      const_cast<volVectorField &>(
          forceFieldPtr_())
          .correctBoundaryConditions();
    }
  }

  if (forceFieldPtr_.valid()) {
    force_ = forceFieldPtr_().boundaryField()[patch().index()];
  }

  // Convert the force field to a traction field
  // Note: this assumes small strains / linear geometry
  const vectorField traction(force_ / patch().magSf());

  // Apply traction
  // The code below comes from tractionDisplacement:updateCoeffs()

  const dictionary &mechanicalProperties =
      db().lookupObject<IOdictionary>("mechanicalProperties");

  const dictionary &thermalProperties =
      db().lookupObject<IOdictionary>("thermalProperties");

  const fvPatchField<scalar> &rho =
      patch().lookupPatchField<volScalarField, scalar>("rho");

  const fvPatchField<scalar> &rhoE =
      patch().lookupPatchField<volScalarField, scalar>("E");

  const fvPatchField<scalar> &nu =
      patch().lookupPatchField<volScalarField, scalar>("nu");

  const scalarField E(rhoE / rho);
  const scalarField mu(E / (2.0 * (1.0 + nu)));
  scalarField       lambda(nu * E / ((1.0 + nu) * (1.0 - 2.0 * nu)));
  scalarField       threeK(E / (1.0 - 2.0 * nu));

  if (mechanicalProperties.get<bool>("planeStress")) {
    lambda = nu * E / ((1.0 + nu) * (1.0 - nu));
    threeK = E / (1.0 - nu);
  }

  const scalarField twoMuLambda(2 * mu + lambda);

  const vectorField n(patch().nf());

  const fvPatchField<symmTensor> &sigmaD =
      patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");

  gradient() =
      (traction / rho + twoMuLambda * fvPatchField<vector>::snGrad() - (n & sigmaD)) / twoMuLambda;

  if (thermalProperties.get<bool>("thermalStress")) {
    const fvPatchField<scalar> &threeKalpha =
        patch().lookupPatchField<volScalarField, scalar>("threeKalpha");

    const fvPatchField<scalar> &T =
        patch().lookupPatchField<volScalarField, scalar>("T");

    gradient() += n * threeKalpha * T / twoMuLambda;
  }

  fixedGradientFvPatchVectorField::updateCoeffs();
}

void solidDisplacementFoamForceFvPatchVectorField::write(Ostream &os) const
{
  // Bug-fix: courtesy of Michael@UW at https://www.cfd-online.com/Forums/
  // openfoam-cc-toolkits-fluid-structure-interaction/221892-solved-paraview
  // -cant-read-solids-files-duplicate-entries-keyword-value.html#post762325
  // fixedGradientFvPatchVectorField::write(os);
  fvPatchVectorField::write(os);

  if (forceFieldPtr_.valid()) {
    os.writeKeyword("forceField")
        << forceFieldPtr_().name() << token::END_STATEMENT << nl;
  } else {
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
