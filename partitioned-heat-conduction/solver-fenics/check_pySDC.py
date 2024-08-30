"""
a helper script to test the installation of pySDC and provide advice if an error is detected
"""
import sys
try:
    from pySDC.core.Problem import ptype
except ModuleNotFoundError:
    print("""
ERROR:

Something is wrong with your installation of pySDC. Please make sure that you have cloned the GitHub repository of pySDC and you are using the correct version. Use the following command to clone pySDC version 5.5.0:

git clone --branch 5.5.0 https://github.com/Parallel-in-Time/pySDC <<put folder here>>

You then have to add the folder into which you cloned pySDC to your PYTHONPATH:

export PYTHONPATH=$PYTHONPATH:<<put folder here>>

If you are running this tutorial in an venv (e.g. using the run script), please make sure that you also correctly set the PYTHONPATH in your venv.

Refer to https://parallel-in-time.org/pySDC/ for installation instructions.
""")
    sys.exit(1)
    
