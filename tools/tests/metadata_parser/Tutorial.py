from dataclasses import dataclass,field
from pathlib import Path
from typing import List,Tuple,Optional
import glob
import yaml
import itertools

from .Case import Case
from .Component import Components

system_test_dir = Path(__file__).parent.parent
@dataclass
class Tutorial:
    """
    Represents a tutorial with various attributes and methods.
    """

    name: str
    path: Path
    url: str
    participants: List[str]
    cases: List[Case]
    case_combinations : List[Tuple[Case]] = field(init=False)

    def __post_init__(self):
        for case in self.cases:
            case.tutorial = self
        # get all case combinations

        def get_all_possible_case_combinations(tutorial: Tutorial):
            case_combinations = []
            cases_dict = {}
            for participant in tutorial.participants:
                cases_dict[participant] = []
            for case in tutorial.cases:
                cases_dict[case.participant].append(case)
            
            for combination in itertools.product(*[cases_dict[participant] for participant in tutorial.participants]):
                case_combinations.append(combination)
            return case_combinations
        
        self.case_combinations = get_all_possible_case_combinations(self)



    def __hash__(self) -> int:
        return hash(self.path)

    def __repr__(self) -> str:
        """
        Returns a string representation of the Tutorial.
        """
        return f"""\n{self.name}:
        Path: {self.path}
        URL: {self.url}
        Participants: {self.participants}
        Cases: {self.cases}
        """

    def get_all_case_combinations(self) -> List[List[Case]]:
        cases_combinations = []
        # first sort the cases into a dict by participant
        
        return cases_combinations


    def get_cases_by_strings(self, case_names: List[List[str]]):
        """
        Retrieves potential cases lists based on their names

        Args:
            case_names: The list of case combinations

        Returns:
            A dictionary of potential cases per participant.
        """
        potential_cases = {}
        for participant in self.participants:
            potential_cases[participant] = []
            for case in self.cases:
                if (case.participant == participant) and (case.component.name in component_names):
                    potential_cases[participant].append(case)
        return potential_cases

    def get_potential_cases(self, component_names):
        """
        Retrieves potential cases based on specified component names.

        Args:
            component_names: The list of component names.

        Returns:
            A dictionary of potential cases per participant.
        """
        potential_cases = {}
        for participant in self.participants:
            potential_cases[participant] = []
            for case in self.cases:
                if (case.participant == participant) and (case.component.name in component_names):
                    potential_cases[participant].append(case)
        return potential_cases

    def can_be_run_with_components(self, component_names):
        """
        Checks if the tutorial can be run with the specified component names.

        Args:
            component_names: The list of component names.

        Returns:
            True if the tutorial can be run with the specified components, False otherwise.
        """
        potential_cases = self.get_potential_cases(component_names)
        can_be_run = True
        for participant in self.participants:
            if len(potential_cases[participant]) == 0:
                can_be_run = False
        return can_be_run

    @classmethod
    def from_yaml(cls, path, available_components):
        """
        Creates a Tutorial instance from a YAML file.

        Args:
            path: The path to the YAML file.
            available_components: The Components instance containing available components.

        Returns:
            An instance of Tutorial.
        """
        with open(path, 'r') as f:
            data = yaml.safe_load(f)
            name = data['name']
            path = data['path']
            url = data['url']
            participants = data.get('participants', [])
            cases_raw = data.get('cases', {})
            cases = []
            for case_name in cases_raw.keys():
                cases.append(Case.from_dict(case_name, cases_raw[case_name], available_components))
            return cls(name, path, url, participants, cases)




class Tutorials(list):
    """
    Represents a collection of tutorials.
    """

    def __iter__(self):
        return iter(self.tutorials)

    def __getitem__(self, index):
        return self.tutorials[index]

    def __setitem__(self, index, value):
        self.tutorials[index] = value

    def __len__(self):
        return len(self.tutorials)

    def __init__(self,tutorials: List[Tutorial]):
        """
        Initializes the Tutorials instance with a base path and a list of tutorials.

        Args:
            path: The path to the folder containing the tutorial folders.
            tutorials: The list of tutorials.
        """
        self.tutorials = tutorials

    def filter_by_components(self, component_names: List[str]) -> List[Tutorial]:
        """
        Filters the tutorials based on the specified component names.

        Args:
            component_names: The list of component names.

        Returns:
            A list of filtered tutorials.
        """
        tutorials_filtered = []
        for tutorial in self.tutorials:
            if tutorial.can_be_run_with_components(component_names):
                tutorials_filtered.append(tutorial)
        return tutorials_filtered

    def get_by_path(self, path_to_search) -> Optional[Tutorial]:
        """
        Retrieves a Tutorial by its relative path.

        Args:
            path_to_search: The path of the Tutorial to search for.

        Returns:
            The Tutorial with the specified path, or None if not found.
        """
        for tutorial in self.tutorials:
            if tutorial.path == path_to_search:
                return tutorial

        return None

    @classmethod
    def from_path(cls, path):
        """
        Read ins all the metadata.yaml files available in path/*/metadata.yaml

        Args:
            path: The path containing the tutorial folders

        """
        yaml_files = glob.glob(f'{path}/*/metadata.yaml')
        tutorials = []
        available_components = Components.from_yaml(system_test_dir / "components.yaml")
        for yaml_path in yaml_files:
            tut = Tutorial.from_yaml(yaml_path, available_components)
            tutorials.append(tut)
        return cls(tutorials)
