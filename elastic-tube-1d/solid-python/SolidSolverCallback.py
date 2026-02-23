from __future__ import print_function

import os
import subprocess
import sys


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    solver = os.path.join(here, "SolidSolver.py")
    default_configuration = os.path.join(here, "..", "precice-config-callback.xml")

    args = sys.argv[1:]
    if len(args) == 0 or args[0].startswith("-"):
        args = [default_configuration] + args

    command = [sys.executable, solver] + args
    subprocess.check_call(command)


if __name__ == "__main__":
    main()
