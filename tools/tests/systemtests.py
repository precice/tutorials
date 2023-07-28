
import argparse
from pathlib import Path
from systemtests.CmdLineArguments import CmdLineArguments
from systemtests.Systemtest import Systemtest
from systemtests.TestSuite import TestSuites
from metadata_parser.metdata import Tutorials,Case




from paths import PRECICE_TUTORIAL_DIR,PRECICE_TESTS_RUN_DIR,PRECICE_TESTS_DIR

parser = argparse.ArgumentParser(description='systemtest')

# Add an argument for the components
parser.add_argument('--suites', type=str, help='Comma-separated test-suites to execute')
parser.add_argument('--build_args', type=str, help='Comma-separated list of arguments provided to the components like openfoam:2102,pythonbindings:latest')
parser.add_argument('--rundir', type=str, help='Directory to run the systemstests in.',nargs='?', const=PRECICE_TESTS_RUN_DIR, default=PRECICE_TESTS_RUN_DIR)

# Parse the command-line arguments
args = parser.parse_args()
systemtests_to_run = []
available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)


build_args = CmdLineArguments.from_args(args.build_args)
run_directory = Path(args.rundir)
if args.suites:
    test_suites_requested = args.suites.split(',')
    available_testsuites = TestSuites.from_yaml(PRECICE_TESTS_DIR/ "tests.yaml",available_tutorials)
    test_suites_to_execute = []
    for test_suite_requested in test_suites_requested:
        test_suite_found = available_testsuites.get_by_name(test_suite_requested)
        if not test_suite_found:
            print(f"Warning: did not find the testsuite with name {test_suite_requested}")
        else:
            test_suites_to_execute.append(test_suite_found)
    if not test_suites_to_execute:
        raise Exception(f"No matching test suites with names {test_suites_requested} found. Use print_test_suites.py to get an overview")
    # now convert the test_suites into systemtests
    for test_suite in test_suites_to_execute:
        for tutorial,case_list in test_suite.cases_of_tutorial.items():
            for cases in case_list:
                systemtests_to_run.append(Systemtest(tutorial,build_args,cases))

    

if not systemtests_to_run:
    raise Exception("Did not find any Systemtests to execute.")


print(f"About to run the following systemtest in the directory {run_directory}: \n{systemtests_to_run}")

for systemtest in systemtests_to_run:
    result = systemtest.run(run_directory)
    if not result.success:
        print(f"Failed to run {result.systemtest}")
    else:
        print(f"Success running {result.systemtest}")
    