// Defines the exported functions for the DLL application.

#include <FECore/sdk.h>
#include "PreciceCallback.h"
#include "DiHuMaterial.h"

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

FECORE_PLUGIN void GetPluginVersion(int &major, int &minor, int &patch)
{
	major = 1;	
	minor = 0;
	patch = 0;
}

FECORE_PLUGIN void PluginInitialize(FECoreKernel &fecore)
{
	FECoreKernel::SetInstance(&fecore);
	// register adapter callback
   	REGISTER_FECORE_CLASS(PreciceCallback, "precice_callback");
	// register custom material
    	REGISTER_FECORE_CLASS(DiHuMaterial, "DiHuMaterial");
	// register custom contraction
	REGISTER_FECORE_CLASS(DiHuContraction, "DiHuContraction");
}
