import glob
import os
import yaml
import itertools

import subprocess
from jinja2 import Environment, FileSystemLoader


def remove_relative_path_prefix(string):
    if string.startswith("./"):
        return string[2:]
    return string


class Case:
    def __init__(self, name,participant,path, run_cmd, clean_cmd,component,tutorial=None) -> None:
        self.name = name
        self.participant = participant
        self.path = path
        self.run_cmd = run_cmd
        self.clean_cmd = clean_cmd
        self.tutorial = tutorial
        self.component= component
        self.cleaned_path= remove_relative_path_prefix(path)

        # Some sanity checks
        self._sanity_checks()

    def _sanity_checks(self):
        if not self.component:
            print(f"FATAL ERROR: Tried to instantiate {self.name} but failed becase no component was found in components.yaml")
            exit(1)

    def render_service_template(self,params_to_use):
        # prepare the dict to give to the template
        render_dict = {
            'params': params_to_use,
            'tutorial': self.tutorial.path,
            'folder': self.cleaned_path,
            'run': self.run_cmd
        }
        jinja_env = Environment(loader=FileSystemLoader('.'))
        template = jinja_env.get_template(self.component.template)
        return template.render(render_dict)
        # select the right template
    
    @classmethod
    def from_dict(cls,name,dict,available_components):
        participant = dict["particpant"]
        path = dict["directory"]
        run_cmd = dict["run"]
        clean_cmd = dict["clean"]
        
        component = available_components.get_by_name(dict["component"])
        return cls(name,participant,path, run_cmd, clean_cmd,component)



    def __repr__(self) -> str:
        return f"{self.name}"


class Participant:
    def __init__(self, name):
        self.name = name

    def __repr__(self) -> str:
        return f"{self.name}"

    def __eq__(self, other):
        if isinstance(other, Participant):
            return self.name == other.name
        return False


class Component:
    def __init__(self,name,template,needed_parameters) -> None:
        self.name=name
        self.template = template
        self.needed_parameters = needed_parameters


    def __eq__(self, other):
        if isinstance(other, Component):
            return self.name == other.name
        return False


    def __repr__(self) -> str:
        return f"{self.name}, {self.needed_parameters}"



class Components:
    def __init__(self,components) -> None:
        self.components = components

    @classmethod
    def from_yaml(cls, path):
        components = []
        params = None
        with open(path, 'r') as f:
            data = yaml.safe_load(f)
            for component_name in data:
                #params = data[component_name]["params"]
                params = Parameters.from_components_yaml(data[component_name])
                template = data[component_name]["docker-service-template"]
                components.append(Component(component_name,template,params))
        
        return cls(components)
    

    def get_by_name(self,name_to_search):
        for component in self.components:
            if component.name == name_to_search:
                return component
    
        return None


class Tutorial:
    def __init__(self,name,path,url,participants,cases):
        self.name = name
        self.path= path 
        self.url = url
        self.participants = participants or []
        self.cases = cases or []
        for case in cases:
            case.tutorial=self

    def __repr__(self) -> str:
        return f"""\n{self.name}:
        Path: {self.path}
        URL: {self.url}
        participants: {self.participants}
        cases: {self.cases}
        """
        

    def get_potential_cases(self,component_names):
        potential_cases={}
        for participant in self.participants:
            potential_cases[participant] = []
            for case in self.cases:
                if (case.participant == participant) and (case.component.name in component_names):
                    # we found a potential match
                    potential_cases[participant].append(case)
        return potential_cases

    def can_be_run_with_components(self,component_names):
        # go over participants and find if we have a case with the component for each participant
        potential_cases = self.get_potential_cases(component_names)
        can_be_run = True
        for participant in self.participants:
            if len(potential_cases[participant]) == 0:
                can_be_run = False
        
        return can_be_run



    
    @classmethod
    def from_yaml(cls, path,available_components):
        with open(path, 'r') as f:
            data = yaml.safe_load(f) 
            name = data['name']
            path = data['path']
            url = data['url']
            participants = data.get('participants', [])
            cases_raw = data.get('cases', {})
            cases = []
            for case_name in cases_raw.keys():
                cases.append(Case.from_dict(case_name,cases_raw[case_name],available_components))
            
            return cls(name,path,url, participants,cases)
        

class Tutorials:

    def __init__(self,base_path,tutorials):
        self.base_path = base_path
        self.tutorials = tutorials


    def filter_by_components(self,component_names):
        tutorials_filtered = []
        for tutorial in self.tutorials:
            if tutorial.can_be_run_with_components(component_names):
                tutorials_filtered.append(tutorial)

        return tutorials_filtered

    @classmethod
    def from_path(cls, path):
        yaml_files = glob.glob(f'{path}/*/metadata.yaml')
        print(f"parsing the following yaml files. {yaml_files}")
        tutorials = []

        available_components = Components.from_yaml("./components.yaml")

        for yaml_path in yaml_files:
            tut = Tutorial.from_yaml(yaml_path,available_components)
            tutorials.append(tut)
        return cls(path,tutorials)



class Parameter:
    def __init__(self,description,key,value_options=None,default=None,is_git_ref=True):
        self.description = description
        self.key = key
        self.default=default
        self.is_git_ref = is_git_ref
        self.required = False if self.default else True


    def __eq__(self, other):
        if isinstance(other, Parameter):
            return self.key == other.key 
        return False

    def __hash__(self):
        return (hash(self.key))

    def __repr__(self):
        return f"{self.key}" 

class InputParameters:
    def __init__(self,parameters):
        self.parameters=parameters
    
    @classmethod
    def from_args(cls,cmd_args):
        if not cmd_args:
            return cls({})
        
        params_provided = cmd_args.split(",")
        params = {}
        for param in params_provided:
            key,value = param.split(":")
            params[key] = value

        return cls(params)

    def __repr__(self):
        return f"{self.parameters}"


    def contains(self,paramter_key):
        return paramter_key in self.parameters.keys()

    def get(self,paramter_key):
        return self.parameters[paramter_key]


class Parameters:
    def __init__(self,parameters):
        self.parameters = parameters

    @classmethod
    def from_components_yaml(cls,data):
        parameters = []
        for param_name in data['params']:
            params = data['params'][param_name]
            description = params.get('description',f"No description provided for {param_name}")
            key=param_name
            is_git_ref = params.get('git-ref',True)
            default = params.get('default',None)
            value_options = params.get('value_options',None)

            parameters.append(Parameter(description,key,value_options,default,is_git_ref))
        
        return cls(parameters)

    def to_list(self):
        return self.parameters
    
    def __repr__(self):
        return f"{self.parameters}"





class Systemtest:
    def __init__(self,tutorial,input_paramters,cases):
        self.tutorial=tutorial
        self.input_paramters = input_paramters
        self.cases = cases
        self.params_to_use = {}
        self._check_if_all_params_are_supplied()
       

    def _check_if_all_params_are_supplied(self):
        # now load the required paramters for each component in the cases
        needed_parameters = set()
        for case in self.cases:
            needed_parameters.update(case.component.needed_parameters.to_list())

        for needed_param in needed_parameters:
            if self.input_paramters.contains(needed_param.key):
                self.params_to_use[needed_param.key] = self.input_paramters.get(needed_param.key)
            else:
                if needed_param.required:
                    print(f"{needed_param} is needed to be given via --params to instantiate {self.tutorial.name}")
                    exit(1)
                else:
                    self.params_to_use[needed_param.key] = needed_param.default


    def _get_docker_services(self):
        rendered_services = {}
        for case in self.cases:
            rendered_services[case.name] = case.render_service_template(self.params_to_use)
        
        return rendered_services
    def _copy_into_directory(self):
        # We copy the whole tutorial over into a runs/tutorial_name folder and run it
        
    def run(self):
        rendered_services = self._get_docker_services()
        render_dict={
            'tutorial': self.tutorial.path,
            'services': rendered_services
        }
        jinja_env = Environment(loader=FileSystemLoader('.'))
        template = jinja_env.get_template("docker-compose.template.yaml")
        docker_compose_file = template.render(render_dict)
        print(f"Now running {self.tutorial.name} with {self.cases}")
        return self.__run_docker_compose(docker_compose_file)

    def __run_docker_compose(self,docker_compose_content):
        # Write the Docker Compose file to disk

        current_dir = os.getcwd()
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
                print(f"Ran {self.tutorial.name} sucessfully.")
                os.remove("docker-compose.yaml")
            else:
                # Print the stdout and stderr
                print(stdout.decode())
                print(stderr.decode())
                print(f"Systemtests {self.tutorial.name} failed with code {exit_code}:")
                print("Docker Compose process failed with exit code:", exit_code)
            
            os.chdir(current_dir)
   
        except OSError as e:
            os.chdir(current_dir)
            print("Error executing docker-compose command:", e)


    
    def __repr__(self):
        return f"{self.tutorial.name} {self.cases}"
