#include <FEBioMech/FEElasticMaterialPoint.h>
#include <FEBioMech/FETransIsoMooneyRivlin.h>
#include <FEBioMech/FEActiveContractionMaterial.h>

/*
 * OpenDiHu Material Point
 * Adds field m_gamma to FEElasticMaterialPoint
 * in order to couple with OpenDiHus FastMonodomainSolver
 */ 
class DiHuMaterialPoint 
	: public FEElasticMaterialPoint
{
public:
    	DiHuMaterialPoint(FEMaterialPointData *mp = nullptr)
    	        : FEElasticMaterialPoint(mp), m_gamma(0) {}

    	double m_gamma; // coupling vatiable
};

/*
 * OpenDiHu Material
 * Behaves the same as FETransIsoMooneyRivlin Material
 * but uses DiHuMaterialPoint instead of FEElasticMaterialPoint
 */
class DiHuMaterial 
	: public FETransIsoMooneyRivlin
{
public:

    	DiHuMaterial(FEModel *pfem) 
    	        : FETransIsoMooneyRivlin(pfem) {}

    	virtual FEMaterialPointData *CreateMaterialPointData() override;

    	DECLARE_FECORE_CLASS();
};

/*
 * OpenDiHu Active Contraction
 * Implements OpenDiHus active contraction model,
 * see https://opendihu.readthedocs.io/en/latest/settings/muscle_contraction_solver.html
 * Use in combination with DiHuMaterial
 */
class DiHuContraction 
	: public FEActiveContractionMaterial
{
public:
    	DiHuContraction(FEModel *pfem)
    	        : FEActiveContractionMaterial(pfem) {}

    	virtual mat3ds ActiveStress(FEMaterialPoint &mp, const vec3d &a0) override;
    	virtual tens4ds ActiveStiffness(FEMaterialPoint &mp, const vec3d &a0) override;

    	double m_pmax; 				// maximum PK2 active stress
	double m_lamOpt; 			// 1.2 constant in OpenDiHu
	bool m_enableForceLengthRelation; 	// whether f(lam/lam_opt) should be multiplied
 
    	DECLARE_FECORE_CLASS();
};
