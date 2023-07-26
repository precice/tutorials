#include <precice/SolverInterface.hpp>
#include <FECore/sdk.h>
#include <FECore/FECallBack.h>
#include <FECore/Callback.h>
#include <FECore/FENode.h>
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEAnalysis.h>
#include <FECore/DumpMemStream.h>

#define PARTICIPANT_NAME "Muscle"
#define ELEMENT_SET "MusclePart"
#define MESH_NAME "MuscleMesh"
#define READ_DATA "Gamma"
#define WRITE_DATA "Geometry"

class PreciceCallback : public FECallBack {
public:
   	PreciceCallback(FEModel *pfem) : FECallBack(pfem, CB_INIT | CB_UPDATE_TIME | CB_MAJOR_ITERS), dmp(*pfem) {}
    	void ReadData(FEModel *fem);
    	void WriteData(FEModel *fem);
    	void Init(FEModel *fem);
    	bool Execute(FEModel &fem, int nreason);
    	std::pair<int, vector<double>> getRelevantMaterialPoints(FEModel *fem, const std::string &elementName);

protected:
    	precice::SolverInterface *precice = NULL;
    	int dimensions; 		// precice dimensions
    	double dt; 			// solver timestep size
    	double precice_dt; 		// precice timestep size 
    	int numberOfVerticies; 		// number of vertices of muscle 
    	int meshID;			// precice ID of the muscle mesh
    	std::vector<int> vertexIDs;	// vertex IDs of the muscle mesh

	// Checkpoints used for implicit coupling
    	const std::string &coric = precice::constants::actionReadIterationCheckpoint();
    	const std::string &cowic = precice::constants::actionWriteIterationCheckpoint();
    	FEAnalysis *checkPointStep;
    	DumpMemStream dmp;
    	double checkpoint_time = 0;
    	FETimeStepController *checkpointTimeStepController = nullptr;
};


