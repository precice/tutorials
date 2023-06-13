import os
import subprocess
from typing import List,Dict
from jinja2 import Environment, FileSystemLoader
from dataclasses import dataclass,field
import shutil

from metadata_parser import Tutorial, Case
from .CmdLineArguments import CmdLineArguments
@dataclass
class Systemtest:
    """
    Represents a system test by specifing the cases and the corresponding Tutorial
    """

    tutorial: Tutorial
    cmd_line_args: CmdLineArguments
    cases: List[Case]
    params_to_use: Dict[str,str] = field(init=False)

    def __post_init__(self):
        self.__init_args_to_use()

    def __init_args_to_use(self):
        """
        Checks if all required parameters for the realisation of the cases are supplied in the cmdline arguments.
        If a parameter is missing and it's required, an exception is raised.
        Otherwise, the default value is used if available.

        In the end it populates the args_to_use dict

        Raises:
            Exception: If a required parameter is missing.
        """
        self.params_to_use={}
        print(self.params_to_use)
        needed_parameters = set()
        for case in self.cases:
            needed_parameters.update(case.component.parameters)

        for needed_param in needed_parameters:
            if self.cmd_line_args.contains(needed_param.key):
                self.params_to_use[needed_param.key] = self.cmd_line_args.get(needed_param.key)
            else:
                if needed_param.required:
                    raise Exception(f"{needed_param} is needed to be given via --params to instantiate the systemtest for {self.tutorial.name}")
                else:
                    self.params_to_use[needed_param.key] = needed_param.default

    def __get_docker_services(self) -> Dict[str,str]:
        """
        Renders the service templates for each case using the parameters to use.

        Returns:
            A dictionary of rendered services per case name.
        """
        def render_service_template(case: Case, params_to_use: Dict[str,str] ) -> str:
            render_dict = {
                'params': params_to_use,
                'tutorial': case.tutorial.path,
                'folder': case.path,
                'run': case.run_cmd
            }
            jinja_env = Environment(loader=FileSystemLoader('.'))
            template = jinja_env.get_template(case.component.template)
            return template.render(render_dict)

        rendered_services = {}
        for case in self.cases:
            rendered_services[case.name] = render_service_template(case,self.params_to_use)
        return rendered_services


    def __get_docker_compose_file(self):
        rendered_services = self.__get_docker_services()
        render_dict = {
            'tutorial': self.tutorial.path,
            'services': rendered_services
        }
        jinja_env = Environment(loader=FileSystemLoader('.'))
        template = jinja_env.get_template("docker-compose.template.yaml")
        return template.render(render_dict)

    def __copy_tutorial_into_directory(self,run_directory: Path):
        """
        Copies the entire tutorial into a folder to prepare for running.
        """
        destination = run_directory / f'{self.tutorial.path}_{self.cases}'
        try:
            path.mkdir(parents=True, exist_ok=False)
        except FileExistsError:
            print("Folder already exists")
        else:
            print("Folder was created")
        shutil.copytree(src, destination)


    def __run_docker_compose(self, docker_compose_content):
        """
        Writes the Docker Compose file to disk, executes docker-compose up, and handles the process output.

        Args:
            docker_compose_content: The content of the Docker Compose file.
        """
        tutorial_path = f"../{self.tutorial.path}"
        os.chdir(tutorial_path)
        with open("docker-compose.yaml", 'w') as file:
            file.write(docker_compose_content)

        try:
            # Execute docker-compose command
            process = subprocess.Popen(['docker-compose', 'up'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            exit_code = process.wait()

            if exit_code == 0:
                print(f"Ran {self.tutorial.name} successfully.")
                os.remove("docker-compose.yaml")
            else:
                # Print the stdout and stderr
                print(stdout.decode())
                print(stderr.decode())
                print(f"System test {self.tutorial.name} failed with code {exit_code}:")
                print("Docker Compose process failed with exit code:", exit_code)

            os.chdir(current_dir)

        except OSError as e:
            os.chdir(current_dir)
            print("Error executing docker-compose command:", e)

    def __repr__(self):
        return f"{self.tutorial.name} {self.cases}"
    
    
    def run(self,run_directory:Path):
        """
        Runs the system test by generating the Docker Compose file, copying everything into a run folder, and executing docker-compose up.
        """
        self.__copy_tutorial_into_directory(run_directory)
       
        print(f"Now running {self.tutorial.name} with {self.cases}")

    
