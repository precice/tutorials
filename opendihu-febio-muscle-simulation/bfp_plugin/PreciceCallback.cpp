#include "PreciceCallback.h"
#include <FECore/log.h>
#include <FECore/FEMaterialPoint.h>
#include <FECore/FETimeStepController.h>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include "DiHuMaterial.h"

// Get number of material points and their initial position
std::pair<int, vector<double>> PreciceCallback::getRelevantMaterialPoints(FEModel *fem, const std::string &elementName) {
    	vector<double> vertexPositions;
    	int counter = 0;
    	FEElementSet *elementSet = fem->GetMesh().FindElementSet(elementName);
    	if (!elementSet) {
    	    	feLogError((elementName + std::string("ElementSet not found")).c_str());
    	    	throw FEException("ElementSet not found");
    	}
    	for (int i = 0; i < elementSet->Elements(); i++) { FEElement &element = elementSet->Element(i); for (int j = 0; j < element.GaussPoints(); j++) {
    	        	// iterate over all materialpoints and add initial position to vector
    	        	FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
    	        	vec3d coord = materialPoint->m_r0; 
    	        	vertexPositions.push_back(coord.x);
    	        	vertexPositions.push_back(coord.y);
    	        	vertexPositions.push_back(coord.z);
    	        	counter++;
    	    	}
    	}
    	return std::pair(counter, vertexPositions);
}

// Initialize the precice adapter
void PreciceCallback::Init(FEModel *fem) {
    	feLogInfo("PreciceCallback::Init");

	// Get config path from envrironment
	const char *config = getenv("BFP_CONFIG");
	if (!config) {
		config = "./precice-config.xml";
	}

    	// initialize precice
    	this->precice = new precice::SolverInterface(PARTICIPANT_NAME, config, 0, 1);
    	this->dimensions = this->precice->getDimensions();
    	
    	// Get material point positions
    	FEMesh &femMesh = fem->GetMesh();
    	std::pair<int, vector<double>> vertexInfo = this->getRelevantMaterialPoints(fem, ELEMENT_SET);
    	this->numberOfVerticies = vertexInfo.first;
    	vector<double> vertexPositions = vertexInfo.second;

    	// Initialize precice mesh
    	this->vertexIDs.resize(this->numberOfVerticies);
    	this->meshID = this->precice->getMeshID(MESH_NAME);
    	this->precice->setMeshVertices(this->meshID, this->numberOfVerticies, vertexPositions.data(), this->vertexIDs.data());

    	// Finish initializing precice
    	this->precice_dt = precice->initialize();
    	this->precice->initializeData();

    	feLogInfo("Finished PreciceCallback::Init");
}

bool PreciceCallback::Execute(FEModel &fem, int nreason) {
    	feLogInfo("PreciceCallback::Execute");

    	if (nreason == CB_INIT) {
    	    	this->Init(&fem);
    	} else if (nreason == CB_UPDATE_TIME) {
    	    	if (this->precice->isActionRequired(this->cowic)) {
    	    	    	feLogInfo("CB_UPDATE_TIME - Saving Checkpoint\n");
    	    	    	// Save
    	    	    	// this uses dmp.open(true,true) which leads to the time controller not beeing serialized
    	    	    	// Also setting dmp.open(true,false) leads to segfault dont know why yet
    	    	    	// Switch time controller
    	    	    	delete this->checkpointTimeStepController;
    	    	    	this->checkpointTimeStepController = fem.GetCurrentStep()->m_timeController;
    	    	    	FETimeStepController *newTimeController = new FETimeStepController(&fem);
    	    	    	newTimeController->SetAnalysis(fem.GetCurrentStep());
    	    	    	newTimeController->CopyFrom(this->checkpointTimeStepController);
    	    	    	fem.GetCurrentStep()->m_timeController = newTimeController;
    	    	    	this->checkpoint_time = fem.GetTime().currentTime;
    	    	    	this->dmp.clear();
    	    	    	fem.Serialize(this->dmp);
    	    	    	this->precice->markActionFulfilled(this->cowic);
    	    	}
    	    	// advance timestep
    	    	this->dt = min(this->precice_dt, fem.GetCurrentStep()->m_dt);
    	    	feLogInfo("Current Simulation Time %f\n", fem.GetTime().currentTime);
    	    	feLogInfo("Timestep %f\n", this->dt);
    	    	fem.GetCurrentStep()->m_dt = this->dt;
    	} else if (nreason == CB_MAJOR_ITERS) {
    	    	if (this->precice->isCouplingOngoing()) {
    	    	    	// Read and write precice data
    	    	    	this->ReadData(&fem);
    	    	    	this->WriteData(&fem);
    	    	    	this->precice_dt = this->precice->advance(this->dt);
    	    	    	if (this->precice->isActionRequired(this->coric)) {
    	    	    	    	feLogInfo("CB_MAJOR_ITERS - Restoring Checkpoint\n");
    	    	    	    	// Restore
    	    	    	    	// taken from FEAnalysis.cpp Line 475 ff
    	    	    	    	// restore the previous state
    	    	    	    	this->dmp.Open(false, true); // This does not restore the time controller only if bshallow is false
    	    	    	    	fem.Serialize(this->dmp);
    	    	    	    	FETimeStepController *newTimeController = new FETimeStepController(&fem);
    	    	    	    	newTimeController->SetAnalysis(fem.GetCurrentStep());
    	    	    	    	newTimeController->CopyFrom(this->checkpointTimeStepController);
    	    	    	    	fem.GetCurrentStep()->m_timeController = newTimeController;
    	    	    	    	fem.GetTime().currentTime = this->checkpoint_time;
    	    	    	    	fem.GetCurrentStep()->m_ntimesteps--; // Decrease number of steps because it gets increased right after this
    	    	    	    	this->precice->markActionFulfilled(this->coric);
    	    	    	}
    	    	}
    	} else if (nreason == CB_SOLVED) {
    	    	this->precice->finalize();
    	    	delete precice;
    	}
    	feLogInfo("Finished PreciceCallback::Execute");
    	return true;
}

// Read Data from precice to febio
void PreciceCallback::ReadData(FEModel *fem) {
    	if (this->precice->isReadDataAvailable()) {
    		feLogInfo("PreciceCallback::ReadData");

    	    	// Read data from precice
    	    	std::vector<double> data(this->numberOfVerticies);
    	    	const int dataID = this->precice->getDataID(READ_DATA, this->meshID); 
    	    	precice->readBlockScalarData(dataID, this->numberOfVerticies, this->vertexIDs.data(), data.data());

    	    	// Write data to febio
    	    	int counter = 0;
    	    	FEElementSet *elementSet = fem->GetMesh().FindElementSet(ELEMENT_SET);
    	    	for (int i = 0; i < elementSet->Elements(); i++) {
    	    	    	FEElement &element = elementSet->Element(i);
    	    	    	for (int j = 0; j < element.GaussPoints(); j++) {
    	    	    	    	FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
				auto elastic = materialPoint->ExtractData<DiHuMaterialPoint>();
				if (elastic == nullptr) {
    	    				throw FEException("MaterialPoint is not an instance of DiHuMaterialPoint");
				}
				elastic->m_gamma = data[counter];
    	    	    	    	counter++;
    	    	    	}
    	    	}
    		feLogInfo("Finished PreciceCallback::ReadData");
    	}
}

// Write data from precice to febio
void PreciceCallback::WriteData(FEModel *fem) {
    if (precice->isWriteDataRequired(dt)) {
    	feLogInfo("PreciceCallback::WriteData");

	// Read data from febio
	std::vector<double> data(this->numberOfVerticies * this->dimensions);
        int counter = 0;
        FEElementSet *elementSet = fem->GetMesh().FindElementSet(ELEMENT_SET);
        for (int i = 0; i < elementSet->Elements(); i++) {
            FEElement &element = elementSet->Element(i);
            for (int j = 0; j < element.GaussPoints(); j++) {
                FEMaterialPoint *materialPoint = element.GetMaterialPoint(j);
		auto rt = materialPoint->m_rt;
		data[counter*dimensions] = rt.x;
		data[counter*dimensions + 1] = rt.y;
		data[counter*dimensions + 2] = rt.z;
                counter++;
            }
        }

	// Write data to precice
	const int dataID = this->precice->getDataID(WRITE_DATA, this->meshID); 
        this->precice->writeBlockVectorData(dataID, this->numberOfVerticies, this->vertexIDs.data(), data.data());
    	feLogInfo("Finished PreciceCallback::WriteData");
    }
}
