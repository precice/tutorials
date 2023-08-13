
import argparse
from metadata_parser.metdata import Tutorials, ReferenceResult
from systemtests.TestSuite import TestSuites
from systemtests.SystemtestArguments import SystemtestArguments
from systemtests.Systemtest import Systemtest
from pathlib import Path
from typing import List, Tuple, Optional, Dict
from paths import PRECICE_TESTS_DIR, PRECICE_TUTORIAL_DIR
import hashlib
from jinja2 import Environment, FileSystemLoader
import tarfile
from datetime import datetime

from paths import PRECICE_TUTORIAL_DIR, PRECICE_TESTS_RUN_DIR, PRECICE_TESTS_DIR, PRECICE_REL_OUTPUT_DIR


def create_tar_gz(source_folder: Path, output_filename: Path):
    with tarfile.open(output_filename, "w:gz") as tar:
        print(source_folder, output_filename)
        tar.add(source_folder, arcname=output_filename.name.replace(".tar.gz", ""))


def render_reference_results_info(
        reference_results: List[ReferenceResult],
        arguments_used: SystemtestArguments,
        time: str):
    def calculate_sha1(file_path: Path):
        buffer_size = 65536
        sha1_hash = hashlib.sha1()
        with open(file_path, "rb") as f:
            while True:
                data = f.read(buffer_size)
                if not data:
                    break
                sha1_hash.update(data)
        return sha1_hash.hexdigest()

    files = []
    for reference_result in reference_results:
        files.append({
            'sha1': calculate_sha1(reference_result.path),
            'time': time,
            'name': reference_result.path.name,
        })

    render_dict = {
        'arguments': arguments_used.arguments,
        'files': files
    }
    jinja_env = Environment(loader=FileSystemLoader(PRECICE_TESTS_DIR))
    template = jinja_env.get_template("reference_results.template")
    return template.render(render_dict)


parser = argparse.ArgumentParser(description='generate reference data')
parser.add_argument('--rundir', type=str, help='Directory to run the systemstests in.',
                    nargs='?', const=PRECICE_TESTS_RUN_DIR, default=PRECICE_TESTS_RUN_DIR)
# Parse the command-line arguments
args = parser.parse_args()

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


print(f"About to run the following tests {systemtests_to_run}")
for systemtest in systemtests_to_run:
    systemtest.run_for_reference_results(run_directory)
    reference_result_per_tutorial[systemtest.tutorial] = []

# Put the tar.gz in there
for systemtest in systemtests_to_run:
    reference_result_folder = systemtest.get_system_test_dir() / PRECICE_REL_OUTPUT_DIR
    reference_result_per_tutorial[systemtest.tutorial].append(systemtest.reference_result)
    create_tar_gz(reference_result_folder, systemtest.reference_result.path)

# write readme
for tutorial in reference_result_per_tutorial.keys():
    with open(tutorial.path / "reference_results.md", 'w') as file:
        ref_results_info = render_reference_results_info(
            reference_result_per_tutorial[tutorial], build_args, current_time_string)
        file.write(ref_results_info)
