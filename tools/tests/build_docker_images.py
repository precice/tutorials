
import argparse
from pathlib import Path
from systemtests.SystemtestArguments import SystemtestArguments
from systemtests.Systemtest import Systemtest, display_systemtestresults_as_table
from systemtests.TestSuite import TestSuites
from metadata_parser.metdata import Tutorials
import logging
import time
from paths import PRECICE_TUTORIAL_DIR, PRECICE_TESTS_RUN_DIR, PRECICE_TESTS_DIR


def main():
    parser = argparse.ArgumentParser(description='build docker images')

    # Add an argument for the components
    parser.add_argument('--suites', type=str,
                        help='Comma-separated test-suites to execute')
    parser.add_argument(
        '--build_args',
        type=str,
        help='Comma-separated list of arguments provided to the components like openfoam:2102,pythonbindings:latest')
    parser.add_argument('--rundir', type=str, help='Directory to run the systemstests in.',
                        nargs='?', const=PRECICE_TESTS_RUN_DIR, default=PRECICE_TESTS_RUN_DIR)

    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO', help='Set the logging level')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Configure logging based on the provided log level
    logging.basicConfig(level=args.log_level, format='%(levelname)s: %(message)s')

    print(f"Using log-level: {args.log_level}")

    systemtests_to_run = []
    available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)

    build_args = SystemtestArguments.from_args(args.build_args)
    run_directory = Path(args.rundir)
    if args.suites:
        test_suites_requested = args.suites.split(',')
        available_testsuites = TestSuites.from_yaml(
            PRECICE_TESTS_DIR / "tests.yaml", available_tutorials)
        test_suites_to_execute = []
        for test_suite_requested in test_suites_requested:
            test_suite_found = available_testsuites.get_by_name(
                test_suite_requested)
            if not test_suite_found:
                logging.error(f"Did not find the testsuite with name {test_suite_requested}")
            else:
                test_suites_to_execute.append(test_suite_found)
        if not test_suites_to_execute:
            raise RuntimeError(
                f"No matching test suites with names {test_suites_requested} found. Use print_test_suites.py to get an overview")
        # now convert the test_suites into systemtests
        for test_suite in test_suites_to_execute:
            tutorials = test_suite.cases_of_tutorial.keys()
            for tutorial in tutorials:
                for case, reference_result in zip(
                        test_suite.cases_of_tutorial[tutorial], test_suite.reference_results[tutorial]):
                    systemtests_to_run.append(
                        Systemtest(tutorial, build_args, case, reference_result))

    if not systemtests_to_run:
        raise RuntimeError("Did not find any Systemtests to execute.")

    logging.info(f"About to build the images for the following systemtests:\n {systemtests_to_run}")

    results = []
    for number, systemtest in enumerate(systemtests_to_run):
        logging.info(f"Started building {systemtest},  {number}/{len(systemtests_to_run)}")
        t = time.perf_counter()
        result = systemtest.run_only_build(run_directory)
        elapsed_time = time.perf_counter() - t
        logging.info(f"Building image for {systemtest} took {elapsed_time} seconds")
        results.append(result)

    build_docker_success = True
    for result in results:
        if not result.success:
            logging.error(f"Failed to run {result.systemtest}")
            build_docker_success = False
        else:
            logging.info(f"Success running {result.systemtest}")

    display_systemtestresults_as_table(results)
    if build_docker_success:
        exit(0)
    else:
        exit(1)


if __name__ == '__main__':
    main()
