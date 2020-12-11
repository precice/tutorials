from enum import Enum

class ProblemType(Enum):
    """
    Enum defines problem type. Details see above.
    """
    DIRICHLET = 1  # Dirichlet problem
    NEUMANN = 2  # Neumann problem


class DomainPart(Enum):
    """
    Enum defines which part of the domain [x_left, x_right] x [y_bottom, y_top] we compute.
    """
    LEFT = 1  # left part of domain in simple interface case
    RIGHT = 2  # right part of domain in simple interface case
    CIRCULAR = 3  # circular part of domain in complex interface case
    RECTANGLE = 4  # domain excluding circular part of complex interface case


class Subcycling(Enum):
    """
    Enum defines which kind of subcycling is used
    """
    NONE = 0  # no subcycling, precice_dt == fenics_dt
    MATCHING = 1  # subcycling, where fenics_dt fits into precice_dt, mod(precice_dt, fenics_dt) == 0
    NONMATCHING = 2  # subcycling, where fenics_dt does not fit into precice_dt, mod(precice_dt, fenics_dt) != 0

    # note: the modulo expressions above should be understood in an exact way (no floating point round off problems. For
    # details, see https://stackoverflow.com/questions/14763722/python-modulo-on-floats)