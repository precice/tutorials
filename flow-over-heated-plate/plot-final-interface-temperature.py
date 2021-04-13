import vtk
from matplotlib import pyplot as plt


def vtk_to_dict(case):
    vtkFileName = "solid-{}/preCICE-output/Fluid-Mesh-Solid.final.vtk".format(case)
    #read the vtk file as an unstructured grid
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(vtkFileName)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    # obtain the data
    data = reader.GetOutput()
    n_data = data.GetPointData().GetNumberOfTuples()

    name = "Temperature"
    data_names = []
    i = 0
    max_i = data.GetPointData().GetNumberOfArrays()
    while i < max_i:
        this_data_name = data.GetPointData().GetArray(i).GetName()
        data_names.append(this_data_name)
        if(this_data_name == name):
            data_id = i
            break
        i+=1
        
    data_dict = {}

    if not data_id:
        raise Exception("For file {} name {} not found. Only the following names are available: {}. Aborting!".format(vtkFileName, name, data_names))
    for i in range(n_data):
        data_dict[data.GetPoint(i)] = data.GetPointData().GetArray(data_id).GetValue(i)
    
    return data_dict


cases = []
cases.append('fenics')
cases.append('openfoam')
cases.append('nutils')

case_labels = {'fenics': 'fenics-openfoam', 'openfoam': 'openfoam-openfoam', 'nutils': 'nutils-openfoam', }
styles = [':','-','--']
colors = ['r','b','g','k']
i = 0

for case in cases:
    case_data = vtk_to_dict(case)
    x, t = [p[0] for p in case_data.keys()], list(case_data.values())
    plt.plot(x, t, colors[i%4]+styles[i%3], label=case_labels[case])
    i += 1

plt.ylabel("Temperature")
plt.xlabel("x-coordinate along coupling interface")
plt.legend()
plt.show()
