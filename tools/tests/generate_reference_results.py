import argparse
from metadata_parser.metdata import Tutorials, ReferenceResult
from systemtests.TestSuite import TestSuites
from systemtests.SystemtestArguments import SystemtestArguments
from systemtests.Systemtest import Systemtest, GLOBAL_TIMEOUT, ITERATIONS_LOGS_DIR
from pathlib import Path
from typing import List
import shutil
from paths import PRECICE_TESTS_DIR, PRECICE_TUTORIAL_DIR
import hashlib
from jinja2 import Environment, FileSystemLoader
import tarfile
import subprocess
from datetime import datetime
import logging

from paths import PRECICE_TUTORIAL_DIR, PRECICE_TESTS_RUN_DIR, PRECICE_TESTS_DIR, PRECICE_REL_OUTPUT_DIR
import time


def create_reference_tar_gz(
    system_test_dir: Path,
    exports_dir: Path,
    output_filename: Path,
    iterations_logs: List[tuple[str, Path]],
) -> None:
    """Archive precice-exports and optional iterations logs as separate top-level tar members."""
    stem = output_filename.name.replace(".tar.gz", "")
    exports_staging = system_test_dir / f".{stem}_reference_exports_staging"
    logs_staging = system_test_dir / f".{stem}_reference_logs_staging"
    for staging in (exports_staging, logs_staging):
        if staging.exists():
            shutil.rmtree(staging)
    shutil.copytree(exports_dir, exports_staging)
    try:
        with tarfile.open(output_filename, "w:gz") as tar:
            tar.add(exports_staging, arcname=stem)
            if iterations_logs:
                for rel, src in iterations_logs:
                    dest = logs_staging / rel
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(src, dest)
                tar.add(
                    logs_staging,
                    arcname=f"{stem}.{ITERATIONS_LOGS_DIR}",
                )
    finally:
        shutil.rmtree(exports_staging, ignore_errors=True)
        shutil.rmtree(logs_staging, ignore_errors=True)


def get_machine_informations():
    def command_is_avail(command: str):
        try:
            rc = subprocess.call(['which', command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            return False

        return rc == 0
    uname_info = "uname not available on the machine the systemtests were executed."
    lscpu_info = "lscpu not available on the machine the systemtests were executed."
    if (command_is_avail("uname")):
        result = subprocess.run(["uname", "-a"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            uname_info = result.stdout

    if (command_is_avail("lscpu") and command_is_avail("grep")):
        result_lscpu = subprocess.run(["lscpu"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        result = subprocess.run(["grep", "-v", "Vulner"], input=result_lscpu.stdout,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            lscpu_info = result.stdout

    return (uname_info, lscpu_info)


def render_reference_results_info(
        reference_results: List[ReferenceResult],
        arguments_used: SystemtestArguments,
        time: str):
    def sha256sum(filename):
        # Implementation from https://stackoverflow.com/a/44873382/2254346,
        # compatible with Python 3.10.
        h = hashlib.sha256()
        mv = memoryview(bytearray(128 * 1024))
        with open(filename, 'rb', buffering=0) as f:
            while n := f.readinto(mv):
                h.update(mv[:n])
        return h.hexdigest()

    files = []
    for reference_result in reference_results:
        files.append({
            'sha256': sha256sum(reference_result.path),
            'time': time,
            'name': reference_result.path.name,
        })
    uname, lscpu = get_machine_informations()
    render_dict = {
        'arguments': arguments_used.arguments,
        'files': files,
        'uname': uname,
        'lscpu': lscpu,
    }

    jinja_env = Environment(loader=FileSystemLoader(PRECICE_TESTS_DIR))
    template = jinja_env.get_template("reference_results.metadata.template")
    return template.render(render_dict)


def main():

    parser = argparse.ArgumentParser(description='Generate reference data for systemtests')
    parser.add_argument('--rundir', type=str, help='Directory to run the systemstests in.',
                        nargs='?', const=PRECICE_TESTS_RUN_DIR, default=PRECICE_TESTS_RUN_DIR)
    parser.add_argument('--suites', type=str,
                        help='Comma-separated test suites to generate reference results for. '
                             'If not specified, all suites are used.')
    parser.add_argument('--log_level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO', help='Set the logging level')
    parser.add_argument(
        '--build_args',
        type=str,
        help='Comma-separated overrides for reference_versions.yaml '
             '(e.g., "TUTORIALS_REF:my-feature-branch")')

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level, format='%(levelname)s: %(message)s')

    print(f"Using log_level: {args.log_level}")

    run_directory = Path(args.rundir)

    available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)

    all_test_suites = TestSuites.from_yaml(PRECICE_TESTS_DIR / "tests.yaml", available_tutorials)

    if args.suites:
        test_suites_requested = []
        for name in args.suites.split(','):
            normalized_name = name.strip()
            if normalized_name and normalized_name not in test_suites_requested:
                test_suites_requested.append(normalized_name)

        if not test_suites_requested:
            parser.error(
                "The --suites option did not contain any valid suite names after parsing. "
                "Use print_test_suites.py to get an overview")

        test_suites = []
        unknown_test_suites = []
        for name in test_suites_requested:
            found = all_test_suites.get_by_name(name)
            if not found:
                unknown_test_suites.append(name)
            else:
                test_suites.append(found)

        if unknown_test_suites:
            parser.error(
                f"Unknown test suite name(s): {', '.join(unknown_test_suites)}. "
                "Use print_test_suites.py to get an overview")

        logging.info(f"Filtering to requested suites: {[s.name for s in test_suites]}")
    else:
        test_suites = all_test_suites
        logging.info("No --suites filter specified, generating reference results for all suites.")

    # Read in parameters from reference_versions.yaml, applying any CLI overrides.
    build_args = SystemtestArguments.from_yaml(PRECICE_TESTS_DIR / "reference_versions.yaml")
    overrides = SystemtestArguments.from_args(args.build_args)
    if overrides.arguments:
        build_args = SystemtestArguments(
            {**build_args.arguments, **overrides.arguments})
        logging.info(f"Applied build argument overrides: {overrides}")
    systemtests_to_run = set()

    for test_suite in test_suites:
        tutorials = test_suite.cases_of_tutorial.keys()
        for tutorial in tutorials:
            max_times = test_suite.max_times.get(tutorial, [])
            mtw_list = test_suite.max_time_windows.get(tutorial, [])
            timeouts = test_suite.timeouts.get(tutorial, [])
            run_befores = test_suite.run_befores.get(tutorial, [])
            run_afters = test_suite.run_afters.get(tutorial, [])
            for i, (case, reference_result) in enumerate(zip(
                    test_suite.cases_of_tutorial[tutorial], test_suite.reference_results[tutorial])):
                max_time = max_times[i] if i < len(max_times) else None
                max_time_windows = mtw_list[i] if i < len(mtw_list) else None
                timeout = timeouts[i] if i < len(timeouts) and timeouts[i] is not None else GLOBAL_TIMEOUT
                run_before = run_befores[i] if i < len(run_befores) else None
                run_after = run_afters[i] if i < len(run_afters) else None
                systemtests_to_run.add(
                    Systemtest(
                        tutorial, build_args, case, reference_result,
                        max_time=max_time, max_time_windows=max_time_windows, timeout=timeout,
                        run_before=run_before, run_after=run_after))

    reference_result_per_tutorial = {}
    current_time_string = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    logging.info(f"About to run the following tests {systemtests_to_run}")
    for number, systemtest in enumerate(systemtests_to_run, start=1):
        logging.info(f"Started running {systemtest},  {number}/{len(systemtests_to_run)}")
        t = time.perf_counter()
        result = systemtest.run_for_reference_results(run_directory)
        elapsed_time = time.perf_counter() - t
        logging.info(f"Running {systemtest} took {elapsed_time:^.1f} seconds")
        if not result.success:
            raise RuntimeError(f"Failed to execute {systemtest}")
        reference_result_per_tutorial[systemtest.tutorial] = []

    # Put the tar.gz in there
    for systemtest in systemtests_to_run:
        reference_result_folder = systemtest.get_system_test_dir() / PRECICE_REL_OUTPUT_DIR
        reference_result_per_tutorial[systemtest.tutorial].append(systemtest.reference_result)
        # create folder if needed
        systemtest.reference_result.path.parent.mkdir(parents=True, exist_ok=True)
        if reference_result_folder.exists():
            collected = systemtest._collect_iterations_logs(systemtest.get_system_test_dir())
            create_reference_tar_gz(
                systemtest.get_system_test_dir(),
                reference_result_folder,
                systemtest.reference_result.path,
                collected,
            )
            if collected:
                logging.info(
                    "Archived %d iterations log(s) inside %s",
                    len(collected),
                    systemtest.reference_result.path.name,
                )
        else:
            raise RuntimeError(
                f"Error executing: \n {systemtest} \n Could not find result folder {reference_result_folder}\n Probably the tutorial did not run through properly. Please check corresponding logs")

    # write readme
    for tutorial in reference_result_per_tutorial.keys():
        reference_results_dir = tutorial.path / "reference-results"
        reference_results_dir.mkdir(parents=True, exist_ok=True)
        with open(reference_results_dir / "reference_results.metadata", 'w') as file:
            ref_results_info = render_reference_results_info(
                reference_result_per_tutorial[tutorial], build_args, current_time_string)
            logging.info(f"Writing results for {tutorial.name}")
            file.write(ref_results_info)
    logging.info(f"Done. Please make sure to manually have a look into the reference results before making a PR.")


if __name__ == '__main__':
    main()
