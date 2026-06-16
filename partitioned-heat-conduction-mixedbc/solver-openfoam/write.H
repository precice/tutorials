// We need to disable the runTime.writeTime if condition for implicit coupling
// as the condition only returns true in the very first coupling iteration
// if (runTime.writeTime())
{
    volVectorField gradT(fvc::grad(T));

    volScalarField gradTx(
        IOobject(
            "gradTx",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        gradT.component(vector::X));

    volScalarField gradTy(
        IOobject(
            "gradTy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        gradT.component(vector::Y));

    volScalarField gradTz(
        IOobject(
            "gradTz",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        gradT.component(vector::Z));

    volVectorField DTgradT(
        IOobject(
            "flux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        DT * gradT);
    volScalarField error_total(
        IOobject(
            "error",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        DT * 0);

    // Print error metrics
    Functions::compute_and_print_errors(mesh, alpha, beta, runTime.value());

    // compute and store the error also in a paraView field
    {
        const Foam::volScalarField* T_(&mesh.lookupObject<volScalarField>("T"));

        // Get the locations of the volume centered mesh vertices
        const vectorField& CellCenters = mesh.C();

        for (int i = 0; i < CellCenters.size(); i++)
        {
            const double coord_x = CellCenters[i].x();
            const double coord_y = CellCenters[i].y();
            const double exact_solution = Functions::get_solution(coord_x, coord_y, alpha, beta, runTime.value());

            error_total.ref()[i] = std::abs((exact_solution - T_->internalField()[i]));
        }

        T.correctBoundaryConditions();
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
                const double exact_solution = Functions::get_solution(coord_x, coord_y, alpha, beta, runTime.value());
                error_total.boundaryFieldRef()[j][i] = std::abs(exact_solution - T_->boundaryField()[j][i]);
            }
        }
    }

    runTime.write();
}
