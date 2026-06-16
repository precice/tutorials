"""
Problem setup for partitioned-heat-conduction/fenicsx tutorial
"""
from dolfinx.io import gmsh as gmshio
from dolfinx.mesh import DiagonalType, create_box
import dolfinx.mesh
from my_enums import DomainPart
import numpy as np
from mpi4py import MPI
import gmsh


def get_geometry(domain_part, communicator):
    gmsh.initialize()
    gmsh.model.add("HeatMesh")

    if domain_part is DomainPart.OUTER:
        # In this section the outer box of the domain is defined.
        # It is a cuboid from (0,0,0) to (2,1,2) and has a cuboid cutout from (1, .25, .5) to (2, .75, 1.5)
        # The coupling interface is the surface from the cutout.

        # Definition of the coupling boundary
        def coupling_bc(x):
            tol = 1E-14
            top = np.logical_and.reduce((
                np.isclose(x[1], 0.75, atol=tol),
                np.logical_or(x[0] >= 1, np.isclose(x[0], 1, tol)),
                np.logical_or(x[2] >= 0.5, np.isclose(x[2], 0.5, tol)),
                np.logical_or(x[2] <= 1.5, np.isclose(x[2], 1.5, tol))
            ))

            bot = np.logical_and.reduce((
                np.isclose(x[1], 0.25, tol),
                np.logical_or(x[0] >= 1, np.isclose(x[0], 1, tol)),
                np.logical_or(x[2] >= 0.5, np.isclose(x[2], 0.5, tol)),
                np.logical_or(x[2] <= 1.5, np.isclose(x[2], 1.5, tol))
            ))

            sideCenter = np.logical_and.reduce((
                np.isclose(x[0], 1, tol),
                x[1] >= 0.25,
                x[1] <= 0.75,
                x[2] >= 0.5,
                x[2] <= 1.5
            ))

            sideLeft = np.logical_and.reduce((
                x[0] >= 1,
                x[1] >= 0.25,
                x[1] <= 0.75,
                np.isclose(x[2], 0.5, tol)
            ))

            sideRight = np.logical_and.reduce((
                x[0] >= 1,
                x[1] >= 0.25,
                x[1] <= 0.75,
                np.isclose(x[2], 1.5, tol)
            ))

            return np.logical_or.reduce((top, bot, sideCenter, sideLeft, sideRight))

        # definition of the boundary excluding the coupling boundary

        def boundary_bc(x):
            tol = 1E-14
            or_part = np.logical_or.reduce((
                np.isclose(x[0], 0, tol),
                np.isclose(x[0], 2, tol),
                np.isclose(x[1], 0, tol),
                np.isclose(x[1], 1, tol),
                np.isclose(x[2], 0, tol),
                np.isclose(x[2], 2, tol)
            ))
            and_part = np.logical_and.reduce((
                np.isclose(x[0], 2, tol),
                x[1] >= 0.25,
                x[1] <= 0.75,
                x[2] >= 0.5,
                x[2] <= 1.5,

            ))
            return np.logical_and(or_part, ~and_part)

        # for creating the mesh, gmsh is used
        # create full cuboid
        gmsh.model.occ.addBox(0, 0, 0, 2, 1, 2, 1)
        gmsh.model.occ.synchronize()
        # create cuboid for cutout
        gmsh.model.occ.addBox(1, 0.25, .5, 1, .5, 1, 2)
        gmsh.model.occ.synchronize()

        # remove smaller cuboid from large cuboid
        _, _ = gmsh.model.occ.cut([(3, 1)], [(3, 2)], 3, True, True)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(3, [3], 4, "foo")
        # generate mesh with max. edge length of 0.15
        gmsh.model.mesh.setOrder(2)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.15)
        gmsh.model.mesh.generate(3)
        # convert from gmsh to a dolfinx representation
        outer_mesh = gmshio.model_to_mesh(gmsh.model, communicator, 0, 3).mesh
        gmsh.finalize()
        return outer_mesh, coupling_bc, boundary_bc
    else:
        # inner part of the domain
        # it is defined by a cuboid from (1, .25, .5) to (2, .75, 1.5)

        # 5 sides of the cuboid define the coupling interface
        def coupling_bc(x):
            tol = 1E-14
            return np.logical_or(
                np.logical_or(np.isclose(x[1], 0.75, tol), np.isclose(x[1], 0.25, tol)),
                np.logical_or(np.isclose(x[0], 1, tol),
                              np.logical_or(np.isclose(x[2], 0.5, tol), np.isclose(x[2], 1.5, tol))))

        def boundary_bc(x):
            tol = 1E-14
            return np.isclose(x[0], 2, tol)

        # create cuboid
        gmsh.model.occ.addBox(1, 0.25, .5, 1, .5, 1, 1)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(3, [1], 2, "foo")
        # generate mesh
        gmsh.model.mesh.setOrder(2)
        gmsh.model.mesh.generate(3)
        # convert from gmsh to dolfinx representation
        inner_mesh = gmshio.model_to_mesh(gmsh.model, communicator, 0, 3).mesh
        gmsh.finalize()
        return inner_mesh, coupling_bc, boundary_bc
