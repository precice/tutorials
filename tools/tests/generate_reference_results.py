
import argparse
from metadata_parser.metdata import Tutorials, ReferenceResult
from systemtests.TestSuite import TestSuites
from systemtests.SystemtestArguments import SystemtestArguments
from systemtests.Systemtest import Systemtest
from pathlib import Path
from typing import List
from paths import PRECICE_TESTS_DIR, PRECICE_TUTORIAL_DIR
import hashlib
from jinja2 import Environment, FileSystemLoader
import tarfile
import subprocess
from datetime import datetime
import logging

from paths import PRECICE_TUTORIAL_DIR, PRECICE_TESTS_RUN_DIR, PRECICE_TESTS_DIR, PRECICE_REL_OUTPUT_DIR
import time


def create_tar_gz(source_folder: Path, output_filename: Path):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_folder, arcname=output_filename.name.replace(".tar.gz", ""))


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
        with open(filename, 'rb', buffering=0) as f:
            return hashlib.file_digest(f, 'sha256').hexdigest()

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
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO', help='Set the logging level')

    args = parser.parse_args()

    logging.basicConfig(level=args.log_level, format='%(levelname)s: %(message)s')

    print(f"Using log-level: {args.log_level}")

    run_directory = Path(args.rundir)

    available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)

    test_suites = TestSuites.from_yaml(PRECICE_TESTS_DIR / "tests.yaml", available_tutorials)

    # Read in parameters
    build_args = SystemtestArguments.from_yaml(PRECICE_TESTS_DIR / "reference_versions.yaml")
    systemtests_to_run = set()

    for test_suite in test_suites:
        tutorials = test_suite.cases_of_tutorial.keys()
        for tutorial in tutorials:
            for case, reference_result in zip(
                    test_suite.cases_of_tutorial[tutorial], test_suite.reference_results[tutorial]):
                systemtests_to_run.add(
                    Systemtest(tutorial, build_args, case, reference_result))

    reference_result_per_tutorial = {}
    current_time_string = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    logging.info(f"About to run the following tests {systemtests_to_run}")
    for number, systemtest in enumerate(systemtests_to_run):
        logging.info(f"Started running {systemtest},  {number}/{len(systemtests_to_run)}")
        t = time.perf_counter()
        result = systemtest.run_for_reference_results(run_directory)
        elapsed_time = time.perf_counter() - t
        logging.info(f"Running {systemtest} took {elapsed_time} seconds")
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
            create_tar_gz(reference_result_folder, systemtest.reference_result.path)
        else:
            raise RuntimeError(
                f"Error executing: \n {systemtest} \n Could not find result folder {reference_result_folder}\n Probably the tutorial did not run through properly. Please check corresponding logs")

    # write readme
    for tutorial in reference_result_per_tutorial.keys():
        with open(tutorial.path / "reference_results.metadata", 'w') as file:
            ref_results_info = render_reference_results_info(
                reference_result_per_tutorial[tutorial], build_args, current_time_string)
            logging.info(f"Writing results for {tutorial.name}")
            file.write(ref_results_info)
    logging.info(f"Done. Please make sure to manually have a look into the reference results before making a PR.")


if __name__ == '__main__':
    main()
