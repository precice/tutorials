// This solver is based on previous work of OpenCFD Ltd. In particular, major code
// parts are part of the laplacianFoam solver of OpenFOAM 2012, which served as a
// basis for this solver.
// -------------------------------------------------------------------------------
// Application
//    heatTransfer
//
// Group
//    grpBasicSolvers
//
// Description
//    Modified version of the Laplace equation solver for a scalar quantity with
//    a non-zero RHS.
//
//    \heading Solver details
//    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
//    equation is given by:
//
//    \f[
//        \ddt{T}  = \div \left( D_T \grad T \right) + F
//    \f]
//
//    Where:
//    \vartable
//        T     | Scalar field which is solved for, e.g. temperature
//        D_T   | Diffusion coefficient
//        F     | The RHS which is defined as: beta - 2 - 2 * alpha
//    \endvartable
//
//    \heading Required fields
//    \plaintable
//        T     | Scalar field which is solved for, e.g. temperature
//    \endplaintable


#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    const double alpha = 3;
    const double beta  = 1.2;
    const double rhs   = beta - 2 - 2 * alpha;

    volScalarField f
     (
         IOobject
         (
             "RHS",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         mesh,
         dimensionedScalar(
         "Tdim",
         dimensionSet(0, 0, -1, 1, 0, 0, 0),
         Foam::scalar(rhs))
     );

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T) - fvm::Su(f,T)
             ==
                fvOptions(T)
            );

            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        #include "write.H"

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
