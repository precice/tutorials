
import argparse
from pathlib import Path
from systemtests.SystemtestArguments import SystemtestArguments
from systemtests.Systemtest import Systemtest, GLOBAL_TIMEOUT, display_systemtestresults_as_table
from systemtests.TestSuite import TestSuites
from metadata_parser.metdata import Tutorials, Case
import logging
import time
import os
import sys
from paths import PRECICE_TUTORIAL_DIR, PRECICE_TESTS_RUN_DIR, PRECICE_TESTS_DIR


class _ConsoleLogFormatter(logging.Formatter):
    """Omit level prefix for INFO/DEBUG; keep it for warnings and errors."""

    def format(self, record: logging.LogRecord) -> str:
        if record.levelno >= logging.WARNING:
            return f"{record.levelname}: {record.getMessage()}"
        return record.getMessage()


def main():
    parser = argparse.ArgumentParser(description='systemtest')

    # Add an argument for the components
    parser.add_argument('--suites', type=str,
                        help='Comma-separated test-suites to execute')
    parser.add_argument(
        '--build_args',
        type=str,
        help='Comma-separated list of arguments provided to the components like openfoam:2102,pythonbindings:latest')
    parser.add_argument('--rundir', type=str, help='Directory to run the systemstests in.',
                        nargs='?', const=PRECICE_TESTS_RUN_DIR, default=PRECICE_TESTS_RUN_DIR)

    parser.add_argument('--log_level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO', help='Set the logging level')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Configure logging
    handler = logging.StreamHandler()
    handler.setFormatter(_ConsoleLogFormatter())
    logging.basicConfig(level=args.log_level, handlers=[handler])

    gh_actions = os.environ.get("GITHUB_ACTIONS", "").lower() == "true"
    # Skip ANSI colors when TERM is unset or "dumb" (minimal terminal, common in CI).
    ansi_colors = sys.stdout.isatty() and os.environ.get("TERM", "") not in {"", "dumb"}

    def _style(text: str, color_code: int | None) -> str:
        if not ansi_colors or color_code is None:
            return text
        return f"\x1b[{color_code}m{text}\x1b[0m"

    def _group_start(title: str) -> None:
        if gh_actions:
            print(f"::group::{title}", flush=True)

    def _group_end() -> None:
        if gh_actions:
            print("::endgroup::", flush=True)

    systemtests_to_run = []
    test_suites_to_execute = []
    available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)

    build_args = SystemtestArguments.from_args(args.build_args)
    run_directory = Path(args.rundir)
    if args.suites:
        test_suites_requested = args.suites.split(',')
        available_testsuites = TestSuites.from_yaml(
            PRECICE_TESTS_DIR / "tests.yaml", available_tutorials)
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
                max_times = test_suite.max_times.get(tutorial, [])
                mtw_list = test_suite.max_time_windows.get(tutorial, [])
                timeouts = test_suite.timeouts.get(tutorial, [])
                for i, (case, reference_result) in enumerate(zip(
                        test_suite.cases_of_tutorial[tutorial], test_suite.reference_results[tutorial])):
                    max_time = max_times[i] if i < len(max_times) else None
                    max_time_windows = mtw_list[i] if i < len(mtw_list) else None
                    timeout = timeouts[i] if i < len(timeouts) and timeouts[i] is not None else GLOBAL_TIMEOUT
                    systemtests_to_run.append(
                        Systemtest(tutorial, build_args, case, reference_result, max_time=max_time, max_time_windows=max_time_windows, timeout=timeout))

    if not systemtests_to_run:
        raise RuntimeError("Did not find any Systemtests to execute.")

    total = len(systemtests_to_run)

    if test_suites_to_execute:
        print("Selected test suite(s):", flush=True)
        print(flush=True)
        for test_suite in test_suites_to_execute:
            print(f"- {test_suite.name}", flush=True)
        print(flush=True)

    print(f"About to run {total} test(s) in the directory {run_directory}:", flush=True)
    print(flush=True)
    for number, systemtest in enumerate(systemtests_to_run, start=1):
        print(f"{number}. {systemtest}", flush=True)
    print(flush=True)
    print(f"Using log_level: {args.log_level}", flush=True)

    results = []
    for number, systemtest in enumerate(systemtests_to_run, start=1):
        print(flush=True)
        started_header = f"[{number}/{total}] Started {systemtest}"
        _group_start(started_header)
        try:
            if not gh_actions:
                logging.info(started_header)
            t = time.perf_counter()
            result = systemtest.run(run_directory)
            elapsed_time = time.perf_counter() - t

            if result.success:
                status_label = _style("✅ PASS", 32)
            else:
                status_label = _style("❌ FAIL", 31)
        finally:
            _group_end()

        print(f"{status_label} Finished in {elapsed_time:.1f}s", flush=True)
        print(flush=True)
        results.append(result)

    print(flush=True)
    system_test_success = all(result.success for result in results)

    display_systemtestresults_as_table(results)
    if system_test_success:
        exit(0)
    else:
        exit(1)


if __name__ == '__main__':
    main()
