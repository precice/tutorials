#include "DiHuMaterial.h"

BEGIN_FECORE_CLASS(DiHuMaterial, FETransIsoMooneyRivlin)
END_FECORE_CLASS();

// Create DiHuMaterialPoint instead of FEElasticMaterialPoint
FEMaterialPointData *DiHuMaterial::CreateMaterialPointData() {
    	auto pt = new DiHuMaterialPoint;
    	if (m_ac) pt->Append(m_ac->CreateMaterialPointData());
    	return pt;
}

BEGIN_FECORE_CLASS(DiHuContraction, FEActiveContractionMaterial)
	ADD_PARAMETER(m_pmax, "pmax");
	ADD_PARAMETER(m_lamOpt, "lam_opt");
	ADD_PARAMETER(m_enableForceLengthRelation, "enable_force_length_relation");
END_FECORE_CLASS();

// Uses OpenDiHus active stress calculation
mat3ds DiHuContraction::ActiveStress(FEMaterialPoint &mp, const vec3d &a0) {
	DiHuMaterialPoint &pt = *mp.ExtractData<DiHuMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// calculate dyad of a: AxA = (a x a)
	mat3ds AxA = dyad(a);

	// Calculate active stress using OpenDiHus formula
	double lamRelative = lamd/m_lamOpt;
	double f = 1.0;
	if (m_enableForceLengthRelation
 	    && 0.6 <= lamRelative && lamRelative <= 1.4) {
		f = -25.0/4.0 * lamRelative*lamRelative + 25.0/2.0 * lamRelative - 5.25;
	}
	double saf = 1.0/lamd * m_pmax * f * pt.m_gamma;

	return AxA*saf;
}

// Do not use active stiffness
tens4ds DiHuContraction::ActiveStiffness(FEMaterialPoint &mp, const vec3d &a0) {
	return tens4ds(0.0);
}
