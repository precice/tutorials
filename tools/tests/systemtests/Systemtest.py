import os
import subprocess
from typing import List, Dict, Optional
from jinja2 import Environment, FileSystemLoader
from dataclasses import dataclass, field
import shutil
from pathlib import Path
from paths import PRECICE_REL_OUTPUT_DIR, PRECICE_TOOLS_DIR, PRECICE_REL_REFERENCE_DIR, PRECICE_TESTS_DIR, PRECICE_TUTORIAL_DIR

from metadata_parser.metdata import Tutorial, CaseCombination, Case, ReferenceResult
from .SystemtestArguments import SystemtestArguments

from datetime import datetime
import tarfile

import unicodedata
import re
import logging


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


@dataclass
class FieldCompareResult:
    exit_code: int
    stdout_data: List[str]
    stderr_data: List[str]
    systemtest: Systemtest


@dataclass
class SystemtestResult:
    success: bool
    stdout_data: List[str]
    stderr_data: List[str]
    systemtest: Systemtest


GLOBAL_TIMEOUT = 360


@dataclass
class Systemtest:
    """
    Represents a system test by specifing the cases and the corresponding Tutorial
    """

    tutorial: Tutorial
    arguments: SystemtestArguments
    case_combination: CaseCombination
    reference_result: ReferenceResult
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
        return hash(f"{self.tutorial,self.arguments,self.case_combination}")

    def __post_init__(self):
        self.__init_args_to_use()
        self.env = {}

    def __init_args_to_use(self):
        """
        Checks if all required parameters for the realisation of the cases are supplied in the cmdline arguments.
        If a parameter is missing and it's required, an exception is raised.
        Otherwise, the default value is used if available.

        In the end it populates the args_to_use dict

        Raises:
            Exception: If a required parameter is missing.
        """
        self.params_to_use = {}
        needed_parameters = set()
        for case in self.case_combination.cases:
            needed_parameters.update(case.component.parameters)

        for needed_param in needed_parameters:
            if self.arguments.contains(needed_param.key):
                self.params_to_use[needed_param.key] = self.arguments.get(
                    needed_param.key)
            else:
                if needed_param.required:
                    raise Exception(
                        f"{needed_param} is needed to be given via --params to instantiate the systemtest for {self.tutorial.name}")
                else:
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
                "-C", repository.resolve(),
                "rev-parse",
                "--abbrev-ref" if abbrev_ref else
                "HEAD"], stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True, check=True, timeout=60)
            current_ref = result.stdout.strip()
            return current_ref
        except Exception as e:
            raise RuntimeError(f"An error occurred while getting the current Git ref: {e}") from e

    def _checkout_ref_in_subfolder(self, repository: Path, subfolder: Path, ref: str):
        try:
            result = subprocess.run([
                "git",
                "-C", repository.resolve(),
                "checkout", ref,
                "--", subfolder.resolve()
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
        current_ref = self._get_git_ref(PRECICE_TUTORIAL_DIR)
        ref_requested = self.params_to_use.get("TUTORIALS_REF")
        if ref_requested:
            logging.debug(f"Checking out tutorials {ref_requested} before copying")
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
        except Exception as e:
            logging.debug(f"tools are already copied: {e} ")

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
                # process.send_signal(9)
                raise KeyboardInterrupt from k
            except Exception as e:
                logging.critical(
                    f"Systemtest {self} had serious issues executin the docker compose command about to kill the docker compose command. Please check the logs! {e}")
                process.kill()
                stdout, stderr = process.communicate()
            stdout_data.extend(stdout.decode().splitlines())
            stderr_data.extend(stderr.decode().splitlines())
            process.poll()
            return FieldCompareResult(process.returncode, stdout_data, stderr_data, self)
        except Exception as e:
            logging.CRITICAL("Error executing docker compose command:", e)
            return FieldCompareResult(1, stdout_data, stderr_data, self)

    def _run_tutorial(self):
        """
        Runs precice couple

        Returns:
            A DockerComposeResult object containing the state.
        """
        docker_compose_content = self.__get_docker_compose_file()
        stdout_data = []
        stderr_data = []

        with open(self.system_test_dir / "docker-compose.tutorial.yaml", 'w') as file:
            file.write(docker_compose_content)
        try:
            # Execute docker-compose command
            process = subprocess.Popen(['docker',
                                        'compose',
                                        '--file',
                                        'docker-compose.tutorial.yaml',
                                        'up',
                                        "--build"],
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
                    f"Systemtest {self} had serious issues executin the docker compose command about to kill the docker compose command. Please check the logs! {e}")
                process.kill()
                stdout, stderr = process.communicate()

            stdout_data.extend(stdout.decode().splitlines())
            stderr_data.extend(stderr.decode().splitlines())

            return DockerComposeResult(process.returncode, stdout_data, stderr_data, self)
        except Exception as e:
            logging.critical(f"Error executing docker compose command: {e}")
            return DockerComposeResult(1, stdout_data, stderr_data, self)

    def __repr__(self):
        return f"{self.tutorial.name} {self.case_combination}"

    def __handle_docker_compose_failure(self, result: DockerComposeResult):
        with open(self.system_test_dir / "stdout.log", 'w') as stdout_file:
            stdout_file.write("\n".join(result.stdout_data))
        with open(self.system_test_dir / "stderr.log", 'w') as stderr_file:
            stderr_file.write("\n".join(result.stderr_data))

        logging.critical("Docker Compose failed, skipping fieldcompare and writing stdout.log and stderr.log")

    def __handle_field_compare_failure(self, result: FieldCompareResult):
        with open(self.system_test_dir / "stdout.log", 'w') as stdout_file:
            stdout_file.write("\n".join(result.stdout_data))
        with open(self.system_test_dir / "stderr.log", 'w') as stderr_file:
            stderr_file.write("\n".join(result.stderr_data))
        logging.error("Fieldcompare failed, writing stdout.log and stderr.log")

    def run(self, run_directory: Path):
        """
        Runs the system test by generating the Docker Compose file, copying everything into a run folder, and executing docker-compose up.
        """
        self.__copy_tutorial_into_directory(run_directory)
        self.__copy_tools(run_directory)
        self.__put_gitignore(run_directory)
        std_out: List[str] = []
        std_err: List[str] = []
        uid, gid = self.__get_uid_gid()
        self.env["UID"] = uid
        self.env["GID"] = gid
        self.__write_env_file()
        docker_compose_result = self._run_tutorial()
        std_out.extend(docker_compose_result.stdout_data)
        std_err.extend(docker_compose_result.stderr_data)
        if docker_compose_result.exit_code != 0:
            self.__handle_docker_compose_failure(docker_compose_result)
            return SystemtestResult(False, std_out, std_err, self)

        fieldcompare_result = self._run_field_compare()
        std_out.extend(fieldcompare_result.stdout_data)
        std_err.extend(fieldcompare_result.stderr_data)
        if fieldcompare_result.exit_code != 0:
            self.__handle_field_compare_failure(fieldcompare_result)
            return SystemtestResult(False, std_out, std_err, self)

        # self.__cleanup()
        return SystemtestResult(True, std_out, std_err, self)

    def run_for_reference_results(self, run_directory: Path):
        """
        Runs the system test by generating the Docker Compose files to generate the reference results
        """
        self.__copy_tutorial_into_directory(run_directory)
        self.__copy_tools(run_directory)
        self.__put_gitignore(run_directory)
        std_out: List[str] = []
        std_err: List[str] = []
        uid, gid = self.__get_uid_gid()
        self.env["UID"] = uid
        self.env["GID"] = gid
        self.__write_env_file()
        docker_compose_result = self._run_tutorial()
        std_out.extend(docker_compose_result.stdout_data)
        std_err.extend(docker_compose_result.stderr_data)
        if docker_compose_result.exit_code == 0:
            return SystemtestResult(True, std_out, std_err, self)
        else:
            self.__handle_docker_compose_failure(docker_compose_result)
            return SystemtestResult(False, std_out, std_err, self)

    def get_system_test_dir(self) -> Path:
        return self.system_test_dir
