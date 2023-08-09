import os
import subprocess
from typing import List, Dict, Tuple
from jinja2 import Environment, FileSystemLoader
from dataclasses import dataclass, field
import shutil
from pathlib import Path

from metadata_parser.metdata import Tutorial, CaseCombination,Case
from .CmdLineArguments import CmdLineArguments

from datetime import datetime


import unicodedata
import re


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


system_test_dir = Path(__file__).parent.parent


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


@dataclass
class Systemtest:
    """
    Represents a system test by specifing the cases and the corresponding Tutorial
    """

    tutorial: Tutorial
    cmd_line_args: CmdLineArguments
    case_combination: CaseCombination
    params_to_use: Dict[str, str] = field(init=False)
    env: Dict[str, str] = field(init=False)

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
            if self.cmd_line_args.contains(needed_param.key):
                self.params_to_use[needed_param.key] = self.cmd_line_args.get(
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
        def render_service_template_per_case(case: Case, params_to_use: Dict[str, str]) -> str:
            render_dict = {
                'run_directory': self.run_directory.resolve(),
                'tutorial_folder': self.tutorial_folder,
                'build_arguments': params_to_use,
                'params': params_to_use,
                'case_folder': case.path,
                'run': case.run_cmd
            }
            jinja_env = Environment(loader=FileSystemLoader(system_test_dir))
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
            'tutorial': self.tutorial.path,
            'services': rendered_services,
        }
        jinja_env = Environment(loader=FileSystemLoader(system_test_dir))
        template = jinja_env.get_template("docker-compose.template.yaml")
        return template.render(render_dict)

    def __get_field_compare_compose_file(self):
        render_dict = {
            'run_directory': self.run_directory.resolve(),
            'tutorial_folder': self.tutorial_folder,
        }
        jinja_env = Environment(loader=FileSystemLoader(system_test_dir))
        template = jinja_env.get_template(
            "docker-compose.field_compare.template.yaml")
        return template.render(render_dict)

    def __copy_tutorial_into_directory(self, run_directory: Path):
        """
        Copies the entire tutorial into a folder to prepare for running.
        """
        current_time_string = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.run_directory = run_directory
        self.tutorial_folder = slugify(f'{self.tutorial.path}_{self.case_combination.cases}_{current_time_string}')
        destination = run_directory / self.tutorial_folder
        src = Path(__file__).parent.parent.parent.parent / self.tutorial.path
        self.system_test_dir = destination
        shutil.copytree(src, destination)

    def __copy_tools(self, run_directory: Path):
        destination = run_directory / "tools"
        src = Path(__file__).parent.parent.parent.parent / "tools"
        try:
            shutil.copytree(src, destination)
        except Exception as e:
            print("tools are already copied: ", e)

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
            # Handle the exception if the 'id -u' or 'id -g' commands fail
            # For example, you can return default values or raise an error.
            print("Error getting group and user id: ", e)

    def __write_env_file(self):
        with open(self.system_test_dir / ".env", "w") as env_file:
            for key, value in self.env.items():
                env_file.write(f"{key}={value}\n")

    def _run_field_compare(self):
        """
        Writes the Docker Compose file to disk, executes docker-compose up, and handles the process output.

        Args:
            docker_compose_content: The content of the Docker Compose file.

        Returns:
            A SystemtestResult object containing the state.
        """
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
                                       cwd=self.system_test_dir)

            # Read the output in real-time
            while True:
                output = process.stdout.readline().decode()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    stdout_data.append(output)
                    print(output, end='')

            # Capture remaining output
            stdout, stderr = process.communicate()
            stdout_data.extend(stdout.decode().splitlines())
            stderr_data.extend(stderr.decode().splitlines())

            exit_code = process.wait()
            return FieldCompareResult(exit_code, stdout_data, stderr_data, self)
        except Exception as e:
            print("Error executing docker compose command:", e)
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
                                       cwd=self.system_test_dir)

            # Read the output in real-time
            while True:
                output = process.stdout.readline().decode()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    stdout_data.append(output)
                    print(output, end='')

            # Capture remaining output
            stdout, stderr = process.communicate()
            stdout_data.extend(stdout.decode().splitlines())
            stderr_data.extend(stderr.decode().splitlines())

            exit_code = process.wait()
            return DockerComposeResult(exit_code, stdout_data, stderr_data, self)
        except Exception as e:
            print("Error executing docker compose command:", e)
            return DockerComposeResult(1, stdout_data, stderr_data, self)

    def __repr__(self):
        return f"{self.tutorial.name} {self.case_combination}"

    def __handle_docker_compose_failure(self, result: DockerComposeResult):
        print("Docker Compose failed, skipping fieldcompare")

    def __handle_field_compare_failure(self, result: FieldCompareResult):
        print("Fieldcompare failed")

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
        if docker_compose_result.exit_code == 1:
            self.__handle_docker_compose_failure(docker_compose_result)
            return SystemtestResult(False, std_out, std_err, self)

        fieldcompare_result = self._run_field_compare()
        std_out.extend(fieldcompare_result.stdout_data)
        std_err.extend(fieldcompare_result.stderr_data)
        if fieldcompare_result.exit_code == 1:
            self.__handle_field_compare_failure(fieldcompare_result)
            return SystemtestResult(False, std_out, std_err, self)

        self.__cleanup()
        return SystemtestResult(True, std_out, std_err, self)
