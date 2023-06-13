
import argparse
from pathlib import Path
from systemtests.CmdLineArguments import CmdLineArguments
from systemtests.Systemtest import Systemtest
from metadata_parser.Tutorial import Tutorials
from metadata_parser.Case import Case


tutorials_path= Path(__file__).parent.parent
systemtests_path = Path(__file__).parent

parser = argparse.ArgumentParser(description='systemtest')

# Add an argument for the components
parser.add_argument('--components', type=str, help='Comma-separated list of components to test')
parser.add_argument('--params', type=str, help='Comma-separated list of arguments provided to the components like openfoam:2102,pythonbindings:latest')

# Parse the command-line arguments
args = parser.parse_args()
tutorials_to_run = None

if args.components:
    # run by components:
    # Extract the components as a list
    components = args.components.split(',')

    # build up the list of tests to run
    available_tutorials = Tutorials.from_path(tutorials_path)

    tutorials_to_run = available_tutorials.filter_by_components(components)


if not tutorials_to_run:
    print("no tutorials to run")
    exit(1)


# generate systemtests from tutorials
systemtests = []
run_directory = systemtests_path / "runs" 
for tutorial in tutorials_to_run:
    # find suitable cases:
    cases_to_run = []
    suitable_cases = tutorial.get_potential_cases(components)
    params = CmdLineArguments.from_args(args.params)
    for participant in tutorial.participants:
        cases_to_run.append(suitable_cases[participant][0])# think of something better to select the cases to run here

    systemtests.append(Systemtest(tutorial,params,cases_to_run))
    print(f"About to run the following systemtests {systemtests}")
    for systest in systemtests: 
        systest.run(run_directory)