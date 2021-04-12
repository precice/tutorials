import vtk

#read the vtk file as an unstructured grid
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("example.vtk")
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

if not data_id:
    raise Exception("Name {} not found. Only the following names are available: {}. Aborting!".format(name, data_names))
for i in range(n_data):
    print(data.GetPoint(i))
    print(data.GetPointData().GetArray(data_id).GetValue(i))

