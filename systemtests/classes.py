import glob
import os
import yaml
import itertools


from jinja2 import Environment, FileSystemLoader


class Case:
    def __init__(self, name,participant,path, run_cmd, clean_cmd,tutorial=None) -> None:
        self.name = name
        self.participant = participant
        self.path = path
        self.run_cmd = run_cmd
        self.clean_cmd = clean_cmd
        self.tutorial = tutorial
    
    @classmethod
    def from_dict(cls,name,dict):
        participant = dict["particpant"]
        path = dict["directory"]
        run_cmd = dict["run"]
        clean_cmd = dict["clean"] 
        return cls(name,participant,path, run_cmd, clean_cmd)



    def __repr__(self) -> str:
        return f"{self.name}"


class Participant:
    def __init__(self, name):
        self.name = name



class Component:
    def __init__(self,name,template,requires) -> None:
        self.name=name
        self.template = template
        self.requires = requires

    def __repr__(self) -> str:
        return f"{self.name}"

class Components:
    def __init__(self,components) -> None:
        self.components = components

    @classmethod
    def from_yaml(cls, path):
        components = []
        with open(path, 'r') as f:
            data = yaml.safe_load(f)
            for component_name in data:
                requires = data[component_name]["requires"]
                template = data[component_name]["docker-service-template"]
                components.append(Component(component_name,template,requires))
        
        return cls(components)
    


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
    
    @classmethod
    def from_yaml(cls, path):
        with open(path, 'r') as f:
            data = yaml.safe_load(f) 
            name = data['name']
            path = data['path']
            url = data['url']
            participants = data.get('participants', [])
            cases_raw = data.get('cases', {})
            cases = []
            for case_name in cases_raw.keys():
                cases.append(Case.from_dict(case_name,cases_raw[case_name]))
            
            return cls(name,path,url, participants,cases)
        

class Tutorials:
    
    def __init__(self,base_path,tutorials):
        self.base_path = base_path
        self.tutorials = tutorials

    @classmethod
    def from_path(cls, path):
        yaml_files = glob.glob(f'{path}/*/metadata.yaml')
        print(f"parsing the following yaml files. {yaml_files}")
        tutorials = []

        for yaml_path in yaml_files:
            tut = Tutorial.from_yaml(yaml_path)
            tutorials.append(tut)
        return cls(path,tutorials)

