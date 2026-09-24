Info << "Reading field T\n"
     << endl;

volScalarField T(
    IOobject(
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

Info << "Creating field gradT\n"
     << endl;

volVectorField gradT(
    IOobject(
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("gradT", dimensionSet(0, -1, 0, 1, 0, 0, 0), vector::zero));

Info << "Creating field maggradT\n"
     << endl;

volScalarField maggradT(
    IOobject(
        "maggradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mag(gradT));

Info << "Reading field U\n"
     << endl;

volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

Info << "Reading transportProperties\n"
     << endl;

IOdictionary transportProperties(
    IOobject(
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE));

Info << "Reading diffusivity DT\n"
     << endl;

dimensionedScalar DT("DT", dimViscosity, transportProperties);

#include "createPhi.H"

#include "createFvOptions.H"
