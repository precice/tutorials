//
// This file is based on the OpenFOAM example application scalarTransportFoam
//
// GitLab link: https://gitlab.com/openfoam/core/openfoam/-/tree/master/applications/solvers/basic/scalarTransportFoam
//
// This modified version doesn't expect a static velocity field U, instead, it repomputes it every time step.
// Therefore it allows the preCICE adapter to provide the velocities.

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  argList::addNote(
      "Dynamic scalar transport equation solver.");

#include "addCheckCaseOptions.H"
#include "createMesh.H"
#include "createTime.H"
#include "setRootCaseLists.H"

  simpleControl simple(mesh);

#include "createFields.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info << "\nCalculating scalar transport\n"
       << endl;

  while (simple.loop()) {
    Info << "Time = " << runTime.timeName() << nl << endl;

    Info << "Recompute phi" << endl;
    fvOptions.correct(U);
    phi = fvc::flux(U);
#include "CourantNo.H"

    while (simple.correctNonOrthogonal()) {
      fvScalarMatrix TEqn(
          fvm::ddt(T) + fvm::div(phi, T) - fvm::laplacian(DT, T) ==
          fvOptions(T));

      TEqn.relax();
      fvOptions.constrain(TEqn);
      TEqn.solve();
      fvOptions.correct(T);
    }

    runTime.write();
  }

  Info << "End\n"
       << endl;

  return 0;
}

// ************************************************************************* //
