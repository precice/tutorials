# BFP Plugin for FEBio 

This plugin enables the coupling of a simple muscle simulation using OpenDiHu and FEBio via preCICE.
The preCICE coupling code is based on [an existing preCICE adapter for FEBio](https://github.com/precice/febio-adapter).

## Installation
The following dependencies are required:
- [preCICE](https://precice.org/installation-overview.html)
- [FEBio](https://febio.org/downloads/)

For the adapter to link against the FEBio SKD successfully you have to perform the following steps:
- Include the SDK when installing FEBio
- Copy *FEBioStudio/lib* to *FEBioStudio/sdk/lib*
- Clone the [FEBio git repository](https://github.com/febiosoftware/FEBio) (branch febio4)
- Copy the folders *FEAMR, FEBio3, FEBioLib, FEBioOpt, FEBioPlot, FEBioRVE, FEBioTest, FEBioXML, FEImgLib, NumCore* and *XML* from the repository to *FEBioStudio/sdk/include*
- Set the path in *build.sh* to your FEBio installation path

You should then be able to compile the plugin with CMake or by running the *build.sh* script.
Then you have to add the following lines to your *FEBioStudio/bin/febio.xml*:
```xml
<import>pathToAdapter/build/lib/libBFPPlugin.so</import>
```
You can also load the plugin in FEBioStudio under *FEBio->Manage FEBio Plugins*.

## Usage
For the PreCICE coupling to work you have to add the following lines to your *model.feb* file 
```xml
<Code>
	<callback name="precice_callback"\>
</Code>
```
In addition you have to change the material to use the custom *DiHuMaterial* and *DiHuContraction* classes
```xml
<Material>
	<material id="material_id" name="material_name" type="DiHuMaterial">
		<density>OpenDiHu density</density>
		<k>1000</k>
		<pressure_model>default</pressure_model>
		<c1>OpenDiHu c1</c1>
		<c2>OpenDiHu c2</c2>
		<c3>0</c3>
		<c4>1</c4>
		<c5>OpenDiHu b</c5>
		<lam_max>1</lam_max>
		<fiber type="vector">
			<vector>Direction of fibers</vector>
		</fiber>
		<active_contraction type="DiHuContraction">
			<pmax>OpenDiHu Pmax</pmax>
			<lam_opt>1.2</lam_opt>
			<enable_force_length_relation>1</enable_force_length_relation>
		</active_contraction>
	</material>
</Material>
```
Note that the following names are hardcoded in the plugin:
- *Muscle*: name of the preCICE participant 
- *MuscleMesh*: name of the preCICE mesh
- *Gamma*, *Geometry*: name of the preCICE fields to couple
- *MusclePart*: name of the FEBio part to couple 

By default the adapter will use *./precice-config.xml* for the preCICE configuration path.
To change this, set the *BFP_CONFIG* environment variable.
You can then run the case with
```bash
BFP_CONFIG="config" mpirun -n 1 FEBioStudio/bin/febio4 model.feb
```

