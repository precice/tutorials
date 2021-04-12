import os


def writeOutputToVTK(time, name, dx, data, datanames):

    if not isinstance(data, list):
        data = list(data)
    if not isinstance(datanames, list):
        datanames = list(datanames)

    n_datasets = data.__len__()
    assert n_datasets == datanames.__len__()

    outpath = os.path.join(os.getcwd(), './output')

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    filename = name + str(time) + ".vtk"
    filepath = os.path.join(outpath, filename)

    i = 0

    f = open(filepath, 'w')

    f.write("# vtk DataFile Version 2.0")
    f.write("\n")
    f.write("\n")

    f.write("ASCII")
    f.write("\n")
    f.write("\n")

    f.write("DATASET UNSTRUCTURED_GRID")
    f.write("\n")
    f.write("\n")

    f.write("POINTS ")
    f.write(str(len(data[i])))
    f.write(" float")
    f.write("\n")
    f.write("\n")

    for k in range(len(data[i])):
        f.write(str("{:.16e}".format(0.0 + k * dx)))
        f.write(" 0.0000000000000000e+00 0.0000000000000000e+00")
        f.write("\n")
    f.write("\n")

    f.write("POINT_DATA ")
    f.write(str(len(data[i])))
    f.write("\n")
    f.write("\n")

    for dataname in datanames:

        if (i == 0):
            f.write("VECTORS ")
        else:
            f.write("SCALARS ")

        f.write(dataname)
        f.write(" float")
        f.write("\n")
        if (i != 0):
            f.write("LOOKUP_TABLE default")
            f.write("\n")
        for element in data[i]:
            f.write(str("{:.16e}".format(element)))
            if (i == 0):
                f.write(" 0.0000000000000000e+00 0.0000000000000000e+00")
            f.write("\n")
        f.write("\n")
        f.write("\n")
        i = i + 1
    f.close()
