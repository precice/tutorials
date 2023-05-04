// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include "properties.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/loadsolution.hh>

#include <map>
#include <list>

int main(int argc, char** argv)
{   
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

    // make some types available
    using QuantityType = Dumux::Precice::QuantityType;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);
    Parameters::print();

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    ////////////////////////
    // Initialize preCICE //
    ////////////////////////

    // Tell preCICE about:
    // - Name of solver
    // - What rank of how many ranks this instance is
    // Configure preCICE. For now the config file is hardcoded.
    std::string preciceConfigFilename = "../precice-config-full-sim.xml";
    if (argc > 2)
        preciceConfigFilename = argv[argc - 1];

    auto &couplingInterface = Dumux::Precice::CouplingAdapter::getInstance();
    
    if (getParam<bool>("Precice.RunWithCoupling") == true)
    {
        couplingInterface.announceSolver("Macro-heat", preciceConfigFilename,
                                     mpiHelper.rank(), mpiHelper.size());
        // verify that dimensions match
        const int preciceDim = couplingInterface.getDimensions();
        const int dim = int(leafGridView.dimension); 
        std::cout << "coupling Dims = " << dim << " , leafgrid dims = " << dim << std::endl;
        if (preciceDim != dim)
            DUNE_THROW(Dune::InvalidStateException, "Dimensions do not match");
    }

    // get mesh coordinates 
    std::string meshName = "macro-mesh";
    std::vector<double> coords;  
    std::vector<int> coupledElementIdxs;
    
    // coordinate loop (created vectors are 1D)
    // these positions of cell centers are later communicated to precice
    std::cout << "Coordinates: " << std::endl;
    for (const auto &element : elements(leafGridView))
    {
        auto fvGeometry = localView(*gridGeometry); 
        fvGeometry.bindElement(element);
        for (const auto &scv : scvs(fvGeometry))
        { 
            coupledElementIdxs.push_back(scv.elementIndex());
            const auto &pos = scv.center();
            for (const auto p : pos)
            {
                coords.push_back(p);
                std::cout <<  p << "  ";
            }
            std::cout << " ;" << std::endl;         
        }
    }
    std::cout << "Number of Coupled Cells:" << coupledElementIdxs.size() << std::endl;

    // initialize preCICE
    double preciceDt;
    auto numberOfElements = coords.size()/couplingInterface.getDimensions();
    if (getParam<bool>("Precice.RunWithCoupling") == true)
    {
        preciceDt = couplingInterface.setMeshAndInitialize(
            meshName, numberOfElements, coords);
        // couples between dumux element indices and preciceIndices; 
        couplingInterface.createIndexMapping(coupledElementIdxs);
    }

    // initialize the coupling data
    std::list<std::string> readDataNames = {"k_00", "k_01", "k_10", "k_11", "porosity"};
    std::map<std::string,size_t> readDataIDs;
    std::string writeDataName = "concentration";
    int temperatureID;
    std::vector<double> temperatures;
    std::vector<double> porosities;
    std::vector<double> k_00;
    std::vector<double> k_01;
    std::vector<double> k_10;
    std::vector<double> k_11;
    std::map<std::string, std::vector<double>> readData {{"k_00", k_00}, {"k_01", k_01}, {"k_10", k_10}, {"k_11", k_11}, {"porosity", porosities}};

    if (getParam<bool>("Precice.RunWithCoupling") == true)
    {
        for (auto iter = readDataNames.begin(); iter != readDataNames.end(); iter++)
            readDataIDs[iter->c_str()] = couplingInterface.announceScalarQuantity(iter->c_str()); //Ids are 0 to 4 in order of readDataNames
        temperatureID = couplingInterface.announceScalarQuantity(writeDataName);
        problem->spatialParams().updatePreciceDataIds();
    }

    // get some time loop parameters
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    double dt;
    if (getParam<bool>("Precice.RunWithCoupling") == true)
        dt = preciceDt;
    else
        dt = getParam<Scalar>("TimeLoop.InitialDt");

    // the solution vector (initialized with zeros) NElements x 2(pressure, temperature)
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());          
    problem->applyInitialSolution(x);  
    auto xOld = x;

    // initialize the coupling data
    for (int solIdx=0; solIdx< numberOfElements; ++solIdx)
    {
        temperatures.push_back(x[solIdx][problem->returnTemperatureIdx()]);
    };

    if (getParam<bool>("Precice.RunWithCoupling") == true)
    {
        couplingInterface.writeQuantityVector(temperatureID, temperatures);
        if (couplingInterface.hasToWriteInitialData()){ //not called in our example
            couplingInterface.writeQuantityVector(temperatureID, temperatures);
            couplingInterface.writeQuantityToOtherSolver(temperatureID, QuantityType::Scalar);
            couplingInterface.announceInitialDataWritten();
        }
        couplingInterface.initializeData();

        // read porosity and conductivity data from other solver
        for (auto iter = readDataNames.begin(); iter != readDataNames.end(); iter++){
            couplingInterface.readQuantityFromOtherSolver(readDataIDs[iter->c_str()], QuantityType::Scalar);
            readData[iter->c_str()] = couplingInterface.getQuantityVector(readDataIDs[iter->c_str()]);
            // we use writeScalarQuantityOnFace to make data accessible (element not face wise) within problem 
            for (const auto &i : coupledElementIdxs)
                couplingInterface.writeScalarQuantityOnFace(readDataIDs[iter->c_str()], i, readData[iter->c_str()][i]);
        }
    }
    
    // the grid variables                           
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); 
    // add model specific output fields
    vtkWriter.addField(problem->getPorosity(), "porosity");
    vtkWriter.addField(problem->getK00(), "k00");
    vtkWriter.addField(problem->getK01(), "k01");
    vtkWriter.addField(problem->getK10(), "k10");
    vtkWriter.addField(problem->getK11(), "k11");
    problem->updateVtkOutput(x);
    vtkWriter.write(0.0);

    // output every vtkOutputInterval time step
    const int vtkOutputInterval = getParam<int>("TimeLoop.OutputInterval");

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    int n = 0; //counts timesteps for the output interval
    std::cout << "Time Loop starts" << std::endl;
    timeLoop->start(); do
    {   
        if (getParam<bool>("Precice.RunWithCoupling") == true){
            if (couplingInterface.isCouplingOngoing()== false)
                break;

            // write checkpoint
            if (couplingInterface.hasToWriteIterationCheckpoint()) 
            {
                xOld = x;
                couplingInterface.announceIterationCheckpointWritten();
            }

            // read porosity and conductivity data from other solver
            for (auto iter = readDataNames.begin(); iter != readDataNames.end(); iter++)
            {
                couplingInterface.readQuantityFromOtherSolver(readDataIDs[iter->c_str()], QuantityType::Scalar);
                readData[iter->c_str()] = couplingInterface.getQuantityVector(readDataIDs[iter->c_str()]);
                for (const auto &i : coupledElementIdxs)
                    couplingInterface.writeScalarQuantityOnFace(readDataIDs[iter->c_str()], i, readData[iter->c_str()][i]);
            }
        }
        std::cout << "Solver starts" << std::endl;

        // linearize & solve
        nonLinearSolver.solve(x, *timeLoop);

        for (int solIdx=0; solIdx< numberOfElements; ++solIdx)
            temperatures[solIdx] = x[solIdx][problem->returnTemperatureIdx()];

        if (getParam<bool>("Precice.RunWithCoupling") == true)
        {
            couplingInterface.writeQuantityVector(temperatureID, temperatures);
            couplingInterface.writeQuantityToOtherSolver(temperatureID, QuantityType::Scalar);      
        }

        // advance precice 
        if (getParam<bool>("Precice.RunWithCoupling") == true)
        {
            preciceDt = couplingInterface.advance(dt);
            dt = std::min(preciceDt, std::min(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()), getParam<Scalar>("TimeLoop.MaxDt")));
            if (preciceDt != dt)
            {
                std::cout << "preciceDt too large. We currently assume fixed timestep size but timesteps no longer correspond: preciceDt = " << preciceDt << " and dt =" << dt << std::endl;
                exit(1);
            }
        }
        else
            dt = std::min(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()), getParam<Scalar>("TimeLoop.MaxDt"));
        
        std::cout << "dt: " << dt << std::endl;
        
        if (getParam<bool>("Precice.RunWithCoupling") == true)
        {
            if (couplingInterface.hasToReadIterationCheckpoint()) 
            {
                // make the new solution the old solution
                x = xOld;
                gridVariables->update(x);
                gridVariables->advanceTimeStep();
                couplingInterface.announceIterationCheckpointRead();
            } 
            else //coupling successful
            {   
                n += 1;
                if (n == vtkOutputInterval)
                {
                    problem->updateVtkOutput(x);
                    vtkWriter.write(timeLoop->time());
                    n = 0;
                }
                gridVariables->advanceTimeStep();
                // advance the time loop to the next step
                timeLoop->advanceTimeStep();
                // report statistics of this time step
                timeLoop->reportTimeStep();
            }
        }
        else
        {
            xOld = x;
            gridVariables->advanceTimeStep();
            // advance the time loop to the next step
            timeLoop->advanceTimeStep();
            // report statistics of this time step
            timeLoop->reportTimeStep();
            
            //output every outputinterval steps
            n += 1;
            if (n == vtkOutputInterval)
            {
                problem->updateVtkOutput(x);
                vtkWriter.write(timeLoop->time());
                n = 0;
            }
        }
        // set new dt as suggested by the newton solver or by precice
        timeLoop->setTimeStepSize(dt);
        
        std::cout << "Time: " << timeLoop->time() << std::endl;

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    if (getParam<bool>("Precice.RunWithCoupling") == true){
        couplingInterface.finalize();
    }
    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
