from enum import Enum


class CouplingScheme(Enum):
    SERIAL_FIRST_DIRICHLET = 1
    SERIAL_FIRST_NEUMANN = 2
    PARALLEL = 3

