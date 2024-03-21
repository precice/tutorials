import numpy as np

# Output, geometric parameters
output_file_name = "out.msh"

x_begin, x_end = 0, 1
y_begin, y_end = -0.25, 0
z_begin, z_end = 0, 0.4

# Number of elements. Add one for number of nodes!
n_x, n_y, n_z = 400, 50, 1

# Coordinates in (i, j, k) space. i = 0, ..., n_x etc
xs = np.linspace(x_begin, x_end, n_x + 1)
ys = np.linspace(y_begin, y_end, n_y + 1)
zs = np.linspace(z_begin, z_end, n_z + 1)

# List of 3D points
nodes = []
# Map from (i,j,k) pos to index
indices = dict()

# Build nodes
for i in range(0, n_x + 1):
    for j in range(0, n_y + 1):
        for k in range(0, n_z + 1):
            nodes.append((xs[i], ys[j], zs[k]))
            # CCX index starts at 1!
            v = len(nodes)
            indices[(i, j, k)] = v + 1

print("** Nodes")
print("*Node, NSET=Nall")
for v, (x, y, z) in enumerate(nodes):
    print("{}, {:10.4f}, {:10.4f}, {:10.4f},".format(v + 2, x, y, z))

# Build elements
print("** Volume elements")
print("* Element, TYPE=C3D8, ELSET=Evolumes")

elems_top_surface = []

elem_id = 1
for i in range(0, n_x):
    for j in range(0, n_y):
        for k in range(0, n_z):
            print("{}, {}, {}, {}, {}, {}, {}, {}, {},".format(elem_id,
                                                               indices[(
                                                                   i, j, k)],
                                                               indices[(
                                                                   i + 1, j, k)],
                                                               indices[(
                                                                   i + 1, j + 1, k)],
                                                               indices[(
                                                                   i, j + 1, k)],
                                                               indices[(
                                                                   i, j, k + 1)],
                                                               indices[(
                                                                   i + 1, j, k + 1)],
                                                               indices[(
                                                                   i + 1, j + 1, k + 1)],
                                                               indices[(i, j + 1, k + 1)],))
            if j == n_y - 1:
                elems_top_surface.append(elem_id)
            elem_id += 1


# Set of border. Adapt freely
print("** Nodes, border with y = y_end")
print("*NSET, NSET=Nyend")
str = ""
for i in range(0, n_x + 1):
    for k in range(0, n_z + 1):
        str = "{}, ".format(indices[i, n_y, k])
        print(str)
print("** Nodes, border with y = y_begin")
print("*NSET, NSET=Nybegin")
str = ""
for i in range(0, n_x + 1):
    for k in range(0, n_z + 1):
        str = "{}, ".format(indices[i, 0, k])
        print(str)


# Upper surface
print("*SURFACE, NAME=Sflux_interface")
for id in elems_top_surface:
    print("{}, S5".format(id))
