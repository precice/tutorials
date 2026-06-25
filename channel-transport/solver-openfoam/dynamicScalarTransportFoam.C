//
// This file is based on the OpenFOAM example application scalarTransportFoam
//
// GitLab link: https://gitlab.com/openfoam/core/openfoam/-/tree/master/applications/solvers/basic/scalarTransportFoam
//
// This modified version doesn't expect a static velocity field U, instead, it repomputes it every time step.
// Therefore it allows the preCICE adapter to provide the velocities.

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::addNote(
        "Dynamic scalar transport equation solver.");

#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createDynamicFvMesh.H"

    simpleControl simple(mesh);

#include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nCalculating scalar transport\n"
         << endl;

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        // Info << "Correct U" << endl;
        // fvOptions.correct(U);
        Info << "Recompute phi" << endl;
        phi = fvc::flux(U);
        Info << "CourantNo" << endl;
#include "CourantNo.H"
        Info << "Update mesh" << endl;
        mesh.update();

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn(
                fvm::ddt(T) + fvm::div(phi, T) - fvm::laplacian(DT, T) == fvOptions(T));

            TEqn.relax();
            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        Info << "Compute gradT" << endl;
        gradT = fvc::grad(T);
        maggradT = mag(gradT);
        runTime.write();
    }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
