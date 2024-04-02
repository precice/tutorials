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
// Some helper functions in order to investigate this case
namespace Functions
{
double get_rhs(const double x, const double y, const double alpha, const double beta)
{
    return beta - 2 - 2 * alpha;
}
double get_solution(const double x, const double y, const double alpha, const double beta, const double time)
{
    return 1 + std::pow(x, 2) + (alpha * std::pow(y, 2)) + (beta * time);
}

void compute_and_print_errors(const Foam::fvMesh& mesh, const double alpha, const double beta, const double time)
{
    double error = 0;
    const Foam::volScalarField* T_(&mesh.lookupObject<volScalarField>("T"));

    double max_error = std::numeric_limits<double>::min();

    // Get the locations of the volume centered mesh vertices
    const vectorField& CellCenters = mesh.C();
    unsigned int numDataLocations = CellCenters.size();
    for (int i = 0; i < CellCenters.size(); i++)
    {
        const double coord_x = CellCenters[i].x();
        const double coord_y = CellCenters[i].y();
        const double exact_solution = Functions::get_solution(coord_x, coord_y, alpha, beta, time);
        const auto result = (exact_solution - T_->internalField()[i]) * (exact_solution - T_->internalField()[i]);
        error += result;
        max_error = std::max(result, max_error);
    }

    Info << "\nError metrics at t = " << time << "s:\n";
    Info << "l2 error (sqrt(sum(err^2)/n)):\t\t" << std::sqrt(error / numDataLocations) << endl;
    Info << "Maximum absolute error (max(err)):\t" << std::sqrt(max_error) << endl;
    Info << "Global absolute error (sum(err^2)):\t" << error << "\n";
    Info << endl;
}
} // namespace Functions

int main(int argc, char* argv[])
{
    argList::addNote(
        "Laplace equation solver for a scalar quantity.");

#include "postProcess.H"

#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"

    simpleControl simple(mesh);

#include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nCalculating temperature distribution\n"
         << endl;

    const double alpha = 3;
    const double beta = 1.2;

    // Initialize the RHS with zero
    volScalarField f(
        IOobject(
            "RHS",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedScalar(
            "Tdim",
            dimensionSet(0, 0, -1, 1, 0, 0, 0),
            Foam::scalar(0)));

    // Now assign the respective values to the RHS
    //(strictly speaking not required here, only for space-dependent RHS functions)
    // Get the locations of the volume centered mesh vertices
    const vectorField& CellCenters = mesh.C();
    for (int i = 0; i < CellCenters.size(); i++)
    {
        const double coord_x = CellCenters[i].x();
        const double coord_y = CellCenters[i].y();
        f.ref()[i] = Functions::get_rhs(coord_x, coord_y, alpha, beta);
    }

    for (int j = 0; j < mesh.boundaryMesh().size() - 1; ++j)
    {
        // Get the face centers of the current patch
        const vectorField faceCenters =
            mesh.boundaryMesh()[j].faceCentres();

        // Assign the (x,y,z) locations to the vertices
        for (int i = 0; i < faceCenters.size(); i++)
        {
            const double coord_x = faceCenters[i].x();
            const double coord_y = faceCenters[i].y();
            f.boundaryFieldRef()[j][i] = Functions::get_rhs(coord_x, coord_y, alpha, beta);
        }
    }

    Functions::compute_and_print_errors(mesh, alpha, beta, runTime.value());

    while (simple.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        // We need to update the coded boundary conditions before solving to account for its time dependency properly
        T.correctBoundaryConditions();
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn(
                fvm::ddt(T) - fvm::laplacian(DT, T) - fvm::Su(f, T)
                == fvOptions(T));

            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

#include "write.H"

        runTime.printExecutionTime(Info);
    }

    Info << "End\n"
         << endl;

    return 0;
}
// ************************************************************************* //
