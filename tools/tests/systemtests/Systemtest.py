import hashlib
import json
import subprocess
from typing import List, Dict, Optional, Tuple
from jinja2 import Environment, FileSystemLoader
from dataclasses import dataclass, field
import shutil
from pathlib import Path
from paths import PRECICE_REL_OUTPUT_DIR, PRECICE_TOOLS_DIR, PRECICE_REL_REFERENCE_DIR, PRECICE_TESTS_DIR, PRECICE_TUTORIAL_DIR

ITERATIONS_LOGS_DIR = "iterations-logs"

from metadata_parser.metdata import Tutorial, CaseCombination, Case, ReferenceResult
from .SystemtestArguments import SystemtestArguments

from datetime import datetime
import tarfile
import time

import unicodedata
import re
import logging
import os


GLOBAL_TIMEOUT = 900
SHORT_TIMEOUT = 10


def slugify(value, allow_unicode=False):
    """
    Taken from https://github.com/django/django/blob/master/django/utils/text.py
    Convert to ASCII if 'allow_unicode' is False. Convert spaces or repeated
    dashes to single dashes. Remove characters that aren't alphanumerics,
    underscores, or hyphens. Convert to lowercase. Also strip leading and
    trailing whitespace, dashes, and underscores.
    """
    value = str(value)
    if allow_unicode:
        value = unicodedata.normalize('NFKC', value)
    else:
        value = unicodedata.normalize('NFKD', value).encode(
            'ascii', 'ignore').decode('ascii')
    value = re.sub(r'[^\w\s-]', '', value.lower())
    return re.sub(r'[-\s]+', '-', value).strip('-_')


class Systemtest:
    pass


@dataclass
class DockerComposeResult:
    exit_code: int
    stdout_data: List[str]
    stderr_data: List[str]
    systemtest: Systemtest
    runtime: float  # in seconds


@dataclass
class FieldCompareResult:
    exit_code: int
    stdout_data: List[str]
    stderr_data: List[str]
    systemtest: Systemtest
    runtime: float  # in seconds


@dataclass
class SystemtestResult:
    success: bool
    stdout_data: List[str]
    stderr_data: List[str]
    systemtest: Systemtest
    build_time: float  # in seconds
    solver_time: float  # in seconds
    fieldcompare_time: float  # in seconds


def display_systemtestresults_as_table(results: List[SystemtestResult]):
    """
    Prints the result in a nice tabluated way to get an easy overview
    """
    def _get_length_of_name(results: List[SystemtestResult]) -> int:
        return max(len(str(result.systemtest)) for result in results)

    max_name_length = _get_length_of_name(results)

    header = f"| {'systemtest':<{max_name_length + 2}} "\
        f"| {'success':^7} "\
        f"| {'building time [s]':^17} "\
        f"| {'solver time [s]':^15} "\
        f"| {'fieldcompare time [s]':^21} |"
    separator_plaintext = "+-" + "-" * (max_name_length + 2) + \
        "-+---------+-------------------+-----------------+-----------------------+"
    separator_markdown = "| --- | --- | --- | --- | --- |"

    print(separator_plaintext)
    print(header)
    print(separator_plaintext)

    if "GITHUB_STEP_SUMMARY" in os.environ:
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f:
            print(header, file=f)
            print(separator_markdown, file=f)

    for result in results:
        row = f"| {str(result.systemtest):<{max_name_length + 2}} "\
            f"| {result.success:^7} "\
            f"| {result.build_time:^17.1f} "\
            f"| {result.solver_time:^15.1f} "\
            f"| {result.fieldcompare_time:^21.1f} |"
        print(row)
        print(separator_plaintext)
        if "GITHUB_STEP_SUMMARY" in os.environ:
            with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f:
                print(row, file=f)

    if "GITHUB_STEP_SUMMARY" in os.environ:
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f:
            print("\n\n", file=f)
            print(
                "In case a test fails, download the archive from the bottom of this page and look into each `stdout.log` and `stderr.log`. The time spent in each step might already give useful hints.",
                file=f)
            print(
                "See the [documentation](https://precice.org/dev-docs-system-tests.html#understanding-what-went-wrong).",
                file=f)


@dataclass
class Systemtest:
    """
    Represents a system test by specifing the cases and the corresponding Tutorial
    """

    tutorial: Tutorial
    arguments: SystemtestArguments
    case_combination: CaseCombination
    reference_result: ReferenceResult
    max_time: Optional[float] = None
    params_to_use: Dict[str, str] = field(init=False)
    env: Dict[str, str] = field(init=False)

    def __eq__(self, other) -> bool:
        if isinstance(other, Systemtest):
            return (
                self.tutorial == other.tutorial) and (
                self.arguments == other.arguments) and (
                self.case_combination == other.case_combination)
        return False

    def __hash__(self) -> int:
        return hash(f"{self.tutorial, self.arguments, self.case_combination}")

    def __post_init__(self):
        self.__init_args_to_use()
        self.env = {}

    def __init_args_to_use(self):
        """
        Forwards the command-line arguments to the params_to_use dictionary, substituting any missing arguments with their defaults.

        Previously, this function was also checking if all required parameters were provided, and was raising exceptions for parameters not provided and not having a default value. This check made adding optional parameters with empty defaults (e.g., the TUTORIALS_PR) complicated, and it was removed.
        """

        # Forward all provided arguments to params_to_use
        provided_arguments = self.arguments.arguments
        self.params_to_use = provided_arguments

        # Find out which parameters are needed
        needed_parameters = set()
        for case in self.case_combination.cases:
            needed_parameters.update(case.component.parameters)

        # Substitute defaults for non-provided, needed arguments
        for needed_param in needed_parameters:
            if not needed_param.key in provided_arguments:
                logging.warning(
                    f"No argument provided for needed parameter {needed_param.key}. Substituting with {needed_param.default}")
                self.params_to_use[needed_param.key] = needed_param.default

    def __get_docker_services(self) -> Dict[str, str]:
        """
        Renders the service templates for each case using the parameters to use.

        Returns:
            A dictionary of rendered services per case name.
        """
        try:
            plaform_requested = self.params_to_use.get("PLATFORM")
        except Exception as exc:
            raise KeyError("Please specify a PLATFORM argument") from exc

        self.dockerfile_context = PRECICE_TESTS_DIR / "dockerfiles" / Path(plaform_requested)
        if not self.dockerfile_context.exists():
            raise ValueError(
                f"The path {self.dockerfile_context.resolve()} resulting from argument PLATFORM={plaform_requested} could not be found in the system")

        def render_service_template_per_case(case: Case, params_to_use: Dict[str, str]) -> str:
            render_dict = {
                'run_directory': self.run_directory.resolve(),
                'tutorial_folder': self.tutorial_folder,
                'build_arguments': params_to_use,
                'params': params_to_use,
                'case_folder': case.path,
                'run': case.run_cmd,
                'dockerfile_context': self.dockerfile_context,
            }
            jinja_env = Environment(loader=FileSystemLoader(PRECICE_TESTS_DIR))
            template = jinja_env.get_template(case.component.template)
            return template.render(render_dict)

        rendered_services = {}
        for case in self.case_combination.cases:
            rendered_services[case.name] = render_service_template_per_case(
                case, self.params_to_use)
        return rendered_services

    def __get_docker_compose_file(self):
        rendered_services = self.__get_docker_services()
        render_dict = {
            'run_directory': self.run_directory.resolve(),
            'tutorial_folder': self.tutorial_folder,
            'tutorial': self.tutorial.path.name,
            'services': rendered_services,
            'build_arguments': self.params_to_use,
            'dockerfile_context': self.dockerfile_context,
            'precice_output_folder': PRECICE_REL_OUTPUT_DIR,
        }
        jinja_env = Environment(loader=FileSystemLoader(PRECICE_TESTS_DIR))
        template = jinja_env.get_template("docker-compose.template.yaml")
        return template.render(render_dict)

    def __get_field_compare_compose_file(self):
        render_dict = {
            'run_directory': self.run_directory.resolve(),
            'tutorial_folder': self.tutorial_folder,
            'precice_output_folder': PRECICE_REL_OUTPUT_DIR,
            'reference_output_folder': PRECICE_REL_REFERENCE_DIR + "/" + self.reference_result.path.name.replace(".tar.gz", ""),
        }
        jinja_env = Environment(loader=FileSystemLoader(PRECICE_TESTS_DIR))
        template = jinja_env.get_template(
            "docker-compose.field_compare.template.yaml")
        return template.render(render_dict)

    def _get_git_ref(self, repository: Path, abbrev_ref=False) -> Optional[str]:
        try:
            result = subprocess.run([
                "git",
                "-C", os.fspath(repository.resolve()),
                "rev-parse",
                "--abbrev-ref" if abbrev_ref else
                "HEAD"], stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True, check=True, timeout=60)
            current_ref = result.stdout.strip()
            return current_ref
        except Exception as e:
            raise RuntimeError(f"An error occurred while getting the current Git ref: {e}") from e

    def _fetch_pr(self, repository: Path, pr: str):
        try:
            result = subprocess.run([
                "git",
                "-C", os.fspath(repository.resolve()),
                "fetch",
                "origin",
                "pull/" + pr + "/head"
            ], check=True, timeout=60)

        except Exception as e:
            raise RuntimeError(f"git command returned code {result.returncode}")

    def _fetch_ref(self, repository: Path, ref: str):
        try:
            result = subprocess.run([
                "git",
                "-C", os.fspath(repository.resolve()),
                "fetch"
            ], check=True, timeout=60)
            if result.returncode != 0:
                raise RuntimeError(f"git command returned code {result.returncode}")

        except Exception as e:
            raise RuntimeError(f"An error occurred while fetching origin '{ref}':  {e}")

    def _checkout_ref_in_subfolder(self, repository: Path, subfolder: Path, ref: str):
        try:
            result = subprocess.run([
                "git",
                "-C", os.fspath(repository.resolve()),
                "checkout", ref,
                "--", os.fspath(subfolder.resolve())
            ], check=True, timeout=60)
            if result.returncode != 0:
                raise RuntimeError(f"git command returned code {result.returncode}")

        except Exception as e:
            raise RuntimeError(f"An error occurred while checking out '{ref}' for folder '{repository}': {e}")

    def __copy_tutorial_into_directory(self, run_directory: Path):
        """
        Checks out the requested tutorial ref and copies the entire tutorial into a folder to prepare for running.
        """
        current_time_string = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.run_directory = run_directory
        pr_requested = self.params_to_use.get("TUTORIALS_PR")
        if pr_requested:
            logging.debug(f"Fetching the PR {pr_requested} HEAD reference")
            self._fetch_pr(PRECICE_TUTORIAL_DIR, pr_requested)
        current_ref = self._get_git_ref(PRECICE_TUTORIAL_DIR)
        ref_requested = self.params_to_use.get("TUTORIALS_REF")
        if ref_requested:
            logging.debug(f"Checking out tutorials {ref_requested} before copying")
            self._fetch_ref(PRECICE_TUTORIAL_DIR, ref_requested)
            self._checkout_ref_in_subfolder(PRECICE_TUTORIAL_DIR, self.tutorial.path, ref_requested)

        self.tutorial_folder = slugify(f'{self.tutorial.path.name}_{self.case_combination.cases}_{current_time_string}')
        destination = run_directory / self.tutorial_folder
        src = self.tutorial.path
        self.system_test_dir = destination
        shutil.copytree(src, destination)

        if ref_requested:
            with open(destination / "tutorials_ref", 'w') as file:
                file.write(ref_requested)
            self._checkout_ref_in_subfolder(PRECICE_TUTORIAL_DIR, self.tutorial.path, current_ref)

    def __copy_tools(self, run_directory: Path):
        destination = run_directory / "tools"
        src = PRECICE_TOOLS_DIR
        try:
            shutil.copytree(src, destination)
        except FileExistsError as e:
            logging.debug(f"Tools directory has already been copied to the workspace - skipping.")
        except Exception as e:
            logging.warning(f"Something went wrong while copying the tools directory to the workspace: {e}")

    def __put_gitignore(self, run_directory: Path):
        # Create the .gitignore file with a single asterisk
        gitignore_file = run_directory / ".gitignore"
        with gitignore_file.open("w") as file:
            file.write("*")

    def __cleanup(self):
        shutil.rmtree(self.run_directory)

    def __get_uid_gid(self):
        try:
            uid = int(subprocess.check_output(["id", "-u"]).strip())
            gid = int(subprocess.check_output(["id", "-g"]).strip())
            return uid, gid
        except Exception as e:
            logging.error("Error getting group and user id: ", e)

    def __write_env_file(self):
        with open(self.system_test_dir / ".env", "w") as env_file:
            for key, value in self.env.items():
                env_file.write(f"{key}={value}\n")

    def __unpack_reference_results(self):
        with tarfile.open(self.reference_result.path) as reference_results_tared:
            # specify which folder to extract to
            reference_results_tared.extractall(self.system_test_dir / PRECICE_REL_REFERENCE_DIR)
        logging.debug(
            f"extracting {self.reference_result.path} into {self.system_test_dir / PRECICE_REL_REFERENCE_DIR}")

    def _run_field_compare(self):
        """
        Writes the Docker Compose file to disk, executes docker-compose up, and handles the process output.

        Args:
            docker_compose_content: The content of the Docker Compose file.

        Returns:
            A SystemtestResult object containing the state.
        """
        logging.debug(f"Running fieldcompare for {self}")
        time_start = time.perf_counter()
        self.__unpack_reference_results()
        docker_compose_content = self.__get_field_compare_compose_file()
        stdout_data = []
        stderr_data = []

        with open(self.system_test_dir / "docker-compose.field_compare.yaml", 'w') as file:
            file.write(docker_compose_content)
        try:
            # Execute docker-compose command
            process = subprocess.Popen(['docker',
                                        'compose',
                                        '--file',
                                        'docker-compose.field_compare.yaml',
                                        'up',
                                        '--exit-code-from',
                                        'field-compare'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       start_new_session=True,
                                       cwd=self.system_test_dir)

            try:
                stdout, stderr = process.communicate(timeout=GLOBAL_TIMEOUT)
            except KeyboardInterrupt as k:
                process.kill()
                raise KeyboardInterrupt from k
            except Exception as e:
                logging.critical(
                    f"Systemtest {self} had serious issues executing the docker compose command about to kill the docker compose command. Please check the logs! {e}")
                process.kill()
                process.communicate(timeout=SHORT_TIMEOUT)
            stdout_data.extend(stdout.decode().splitlines())
            stderr_data.extend(stderr.decode().splitlines())
            process.poll()
            elapsed_time = time.perf_counter() - time_start
            return FieldCompareResult(process.returncode, stdout_data, stderr_data, self, elapsed_time)
        except Exception as e:
            logging.CRITICAL("Error executing docker compose command:", e)
            elapsed_time = time.perf_counter() - time_start
            return FieldCompareResult(1, stdout_data, stderr_data, self, elapsed_time)

    @staticmethod
    def _sha256_file(path: Path) -> str:
        """Compute SHA-256 hex digest of a file."""
        h = hashlib.sha256()
        mv = memoryview(bytearray(128 * 1024))
        with open(path, 'rb', buffering=0) as f:
            while n := f.readinto(mv):
                h.update(mv[:n])
        return h.hexdigest()

    def _collect_iterations_logs(
        self, system_test_dir: Path
    ) -> List[Tuple[str, Path]]:
        """
        Collect precice-*-iterations.log files from case dirs.
        Returns list of (relative_path, absolute_path) e.g. ("solid-fenics/precice-Solid-iterations.log", path).
        """
        collected = []
        for case in self.case_combination.cases:
            case_dir = system_test_dir / Path(case.path).name
            if not case_dir.exists():
                continue
            for log_file in case_dir.glob("precice-*-iterations.log"):
                if log_file.is_file():
                    rel = f"{Path(case.path).name}/{log_file.name}"
                    collected.append((rel, log_file))
        return collected

    def __archive_iterations_logs(self):
        """
        Copy precice-*-iterations.log from case dirs into iterations-logs/
        so they are available in CI artifacts (issue #440).
        """
        collected = self._collect_iterations_logs(self.system_test_dir)
        if not collected:
            return
        dest_dir = self.system_test_dir / ITERATIONS_LOGS_DIR
        dest_dir.mkdir(exist_ok=True)
        for rel, src in collected:
            dest_name = Path(rel).name
            if len(collected) > 1:
                prefix = Path(rel).parent.name + "_"
                dest_name = prefix + dest_name
            shutil.copy2(src, dest_dir / dest_name)
        logging.debug(f"Archived {len(collected)} iterations log(s) to {dest_dir} for {self}")

    def __compare_iterations_hashes(self) -> bool:
        """
        Compare current iterations.log hashes against reference sidecar.
        Returns True if comparison passes (or is skipped). Returns False if hashes differ.
        """
        sidecar = self.reference_result.path.with_suffix(".iterations-hashes.json")
        if not sidecar.exists():
            return True
        try:
            ref_hashes = json.loads(sidecar.read_text())
        except (json.JSONDecodeError, OSError) as e:
            logging.warning(f"Could not read iterations hashes from {sidecar}: {e}")
            return True
        if not ref_hashes:
            return True
        collected = self._collect_iterations_logs(self.system_test_dir)
        current = {rel: self._sha256_file(p) for rel, p in collected}
        for rel, expected in ref_hashes.items():
            if rel not in current:
                logging.critical(
                    f"Missing iterations log {rel} (expected from reference); {self} fails"
                )
                return False
            if current[rel] != expected:
                logging.critical(
                    f"Hash mismatch for {rel} (iterations.log regression); {self} fails"
                )
                return False
        if len(current) != len(ref_hashes):
            extra = set(current) - set(ref_hashes)
            logging.critical(
                f"Unexpected iterations log(s) {extra}; {self} fails"
            )
            return False
        return True

    def _build_docker(self):
        """
        Builds the docker image
        """
        logging.debug(f"Building docker image for {self}")
        time_start = time.perf_counter()
        docker_compose_content = self.__get_docker_compose_file()
        with open(self.system_test_dir / "docker-compose.tutorial.yaml", 'w') as file:
            file.write(docker_compose_content)

        stdout_data = []
        stderr_data = []

        try:
            # Execute docker-compose command
            process = subprocess.Popen(['docker',
                                        'compose',
                                        '--file',
                                        'docker-compose.tutorial.yaml',
                                        'build'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       start_new_session=True,
                                       cwd=self.system_test_dir)

            try:
                stdout, stderr = process.communicate(timeout=GLOBAL_TIMEOUT)
            except KeyboardInterrupt as k:
                process.kill()
                # process.send_signal(9)
                raise KeyboardInterrupt from k
            except Exception as e:
                logging.critical(
                    f"systemtest {self} had serious issues building the docker images via the `docker compose build` command. About to kill the docker compose command. Please check the logs! {e}")
                process.communicate(timeout=SHORT_TIMEOUT)
                process.kill()

            stdout_data.extend(stdout.decode().splitlines())
            stderr_data.extend(stderr.decode().splitlines())
            elapsed_time = time.perf_counter() - time_start
            return DockerComposeResult(process.returncode, stdout_data, stderr_data, self, elapsed_time)
        except Exception as e:
            logging.critical(f"Error executing docker compose build command: {e}")
            elapsed_time = time.perf_counter() - time_start
            return DockerComposeResult(1, stdout_data, stderr_data, self, elapsed_time)

    def _run_tutorial(self):
        """
        Runs precice couple

        Returns:
            A DockerComposeResult object containing the state.
        """
        logging.debug(f"Running tutorial {self}")
        time_start = time.perf_counter()
        stdout_data = []
        stderr_data = []
        try:
            # Execute docker-compose command
            process = subprocess.Popen(['docker',
                                        'compose',
                                        '--file',
                                        'docker-compose.tutorial.yaml',
                                        'up'],
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE,
                                       start_new_session=True,
                                       cwd=self.system_test_dir)

            try:
                stdout, stderr = process.communicate(timeout=GLOBAL_TIMEOUT)
            except KeyboardInterrupt as k:
                process.kill()
                # process.send_signal(9)
                raise KeyboardInterrupt from k
            except Exception as e:
                logging.critical(
                    f"Systemtest {self} had serious issues executing the docker compose command about to kill the docker compose command. Please check the logs! {e}")
                process.kill()
                stdout, stderr = process.communicate(timeout=SHORT_TIMEOUT)
                process.kill()

            stdout_data.extend(stdout.decode().splitlines())
            stderr_data.extend(stderr.decode().splitlines())
            elapsed_time = time.perf_counter() - time_start
            return DockerComposeResult(process.returncode, stdout_data, stderr_data, self, elapsed_time)
        except Exception as e:
            logging.critical(f"Error executing docker compose up command: {e}")
            elapsed_time = time.perf_counter() - time_start
            return DockerComposeResult(1, stdout_data, stderr_data, self, elapsed_time)

    def __repr__(self):
        return f"{self.tutorial.name} {self.case_combination}"

    def __write_logs(self, stdout_data: List[str], stderr_data: List[str]):
        with open(self.system_test_dir / "stdout.log", 'w') as stdout_file:
            stdout_file.write("\n".join(stdout_data))
        with open(self.system_test_dir / "stderr.log", 'w') as stderr_file:
            stderr_file.write("\n".join(stderr_data))

    def __apply_precice_max_time_override(self):
        """
        If max_time is set, override <max-time value="..."> in precice-config.xml
        of the copied tutorial directory. Applies to both test runs and reference generation.
        Uses a precise pattern to target only <max-time> tags.
        """
        if self.max_time is None:
            return
        if not (isinstance(self.max_time, (int, float)) and self.max_time > 0):
            logging.warning(
                f"Invalid max_time {self.max_time} for {self}; must be a positive number. Skipping override.")
            return
        config_path = self.system_test_dir / "precice-config.xml"
        if not config_path.exists():
            logging.warning(
                f"Requested max_time override for {self}, but no precice-config.xml "
                f"found in {self.system_test_dir}"
            )
            return
        try:
            text = config_path.read_text()
        except Exception as e:
            logging.warning(f"Could not read {config_path} to apply max_time override: {e}")
            return
        # Target only <max-time value="..."> to avoid modifying time-window-size etc.
        pattern = r'(<max-time\s+value=")([^"]*)(")'
        matches = re.findall(pattern, text)
        if not matches:
            logging.warning(
                f"Requested max_time override for {self}, but no <max-time> tag "
                f"found in {config_path}"
            )
            return
        if len(matches) > 1:
            logging.warning(
                f"Multiple <max-time> tags found in {config_path}; overriding all to {self.max_time}"
            )
        new_text = re.sub(pattern, rf"\g<1>{self.max_time}\g<3>", text)
        try:
            config_path.write_text(new_text)
            logging.info(f"Overwrote max-time in {config_path} to {self.max_time} for {self}")
        except Exception as e:
            logging.warning(f"Failed to write updated {config_path}: {e}")

    def __prepare_for_run(self, run_directory: Path):
        """
        Prepares the run_directory with folders and datastructures needed for every systemtest execution
        """
        self.__copy_tutorial_into_directory(run_directory)
        self.__apply_precice_max_time_override()
        self.__copy_tools(run_directory)
        self.__put_gitignore(run_directory)
        host_uid, host_gid = self.__get_uid_gid()
        self.params_to_use['PRECICE_UID'] = host_uid
        self.params_to_use['PRECICE_GID'] = host_gid

    def run(self, run_directory: Path):
        """
        Runs the system test by generating the Docker Compose file, copying everything into a run folder, and executing docker-compose up.
        """
        self.__prepare_for_run(run_directory)
        std_out: List[str] = []
        std_err: List[str] = []

        docker_build_result = self._build_docker()
        std_out.extend(docker_build_result.stdout_data)
        std_err.extend(docker_build_result.stderr_data)
        if docker_build_result.exit_code != 0:
            self.__write_logs(std_out, std_err)
            logging.critical(f"Could not build the docker images, {self} failed")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=0,
                fieldcompare_time=0)

        docker_run_result = self._run_tutorial()
        std_out.extend(docker_run_result.stdout_data)
        std_err.extend(docker_run_result.stderr_data)
        if docker_run_result.exit_code != 0:
            self.__write_logs(std_out, std_err)
            logging.critical(f"Could not run the tutorial, {self} failed")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=docker_run_result.runtime,
                fieldcompare_time=0)

        self.__archive_iterations_logs()
        if not self.__compare_iterations_hashes():
            self.__write_logs(std_out, std_err)
            logging.critical(
                f"Iterations.log hash comparison failed (regression), {self} failed"
            )
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=docker_run_result.runtime,
                fieldcompare_time=0)

        fieldcompare_result = self._run_field_compare()
        std_out.extend(fieldcompare_result.stdout_data)
        std_err.extend(fieldcompare_result.stderr_data)
        if fieldcompare_result.exit_code != 0:
            self.__write_logs(std_out, std_err)
            logging.critical(f"Fieldcompare returned non zero exit code, therefore {self} failed")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=docker_run_result.runtime,
                fieldcompare_time=fieldcompare_result.runtime)

        # self.__cleanup()
        self.__write_logs(std_out, std_err)
        return SystemtestResult(
            True,
            std_out,
            std_err,
            self,
            build_time=docker_build_result.runtime,
            solver_time=docker_run_result.runtime,
            fieldcompare_time=fieldcompare_result.runtime)

    def run_for_reference_results(self, run_directory: Path):
        """
        Runs the system test by generating the Docker Compose files to generate the reference results
        """
        self.__prepare_for_run(run_directory)
        std_out: List[str] = []
        std_err: List[str] = []
        docker_build_result = self._build_docker()
        std_out.extend(docker_build_result.stdout_data)
        std_err.extend(docker_build_result.stderr_data)
        if docker_build_result.exit_code != 0:
            self.__write_logs(std_out, std_err)
            logging.critical(f"Could not build the docker images, {self} failed")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=0,
                fieldcompare_time=0)

        docker_run_result = self._run_tutorial()
        std_out.extend(docker_run_result.stdout_data)
        std_err.extend(docker_run_result.stderr_data)
        if docker_run_result.exit_code != 0:
            self.__write_logs(std_out, std_err)
            logging.critical(f"Could not run the tutorial, {self} failed")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=docker_run_result.runtime,
                fieldcompare_time=0)

        self.__write_logs(std_out, std_err)
        return SystemtestResult(
            True,
            std_out,
            std_err,
            self,
            build_time=docker_build_result.runtime,
            solver_time=docker_run_result.runtime,
            fieldcompare_time=0)

    def run_only_build(self, run_directory: Path):
        """
        Runs only the build commmand, for example to preheat the caches of the docker builder.
        """
        self.__prepare_for_run(run_directory)
        std_out: List[str] = []
        std_err: List[str] = []
        docker_build_result = self._build_docker()
        std_out.extend(docker_build_result.stdout_data)
        std_err.extend(docker_build_result.stderr_data)
        if docker_build_result.exit_code != 0:
            self.__write_logs(std_out, std_err)
            logging.critical(f"Could not build the docker images, {self} failed")
            return SystemtestResult(
                False,
                std_out,
                std_err,
                self,
                build_time=docker_build_result.runtime,
                solver_time=0,
                fieldcompare_time=0)

        self.__write_logs(std_out, std_err)
        return SystemtestResult(
            True,
            std_out,
            std_err,
            self,
            build_time=docker_build_result.runtime,
            solver_time=0,
            fieldcompare_time=0)

    def get_system_test_dir(self) -> Path:
        return self.system_test_dir
