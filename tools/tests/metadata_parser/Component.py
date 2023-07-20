from dataclasses import dataclass
from typing import Optional,List
import yaml

@dataclass
class Parameter:
    """Represents a parameter needed for the instanciation of a component"""

    description: str
    """The description of the parameter."""
    
    key: str
    """The name of the parameter."""
    
    value_options: Optional[list] = None
    """The optinal list of value options for the parameter. If none is suplied all values are accepted"""
    
    default: Optional[str] = None
    """The default value for the parameter."""
    
    @property
    def required(self) -> bool:
        """
        Check if the parameter need to be supplied via CommandLineArgs

        Returns:
            bool: True if the parameter is required, False otherwise.
        """
        return False if self.default else True
    
    def __eq__(self, other) -> bool:
        if isinstance(other, Parameter):
            return self.key == other.key
        return False

    def __hash__(self) -> int:
        return hash(self.key)

    def __repr__(self) -> str:
        return f"{self.key}"



class Parameters:
    """Represents a collection of parameters."""

    def __init__(self, parameters: List[Parameter]):
        self.parameters  = parameters

    @classmethod
    def from_components_yaml(cls, data):
        """
        Create a list of Paramters from the components YAML data.

        Args:
            data: The components YAML data.
        """
        parameters = []
        for param_name, params in data['build-arguments'].items():
            ## TODO maybe **params
            description = params.get('description', f"No description provided for {param_name}")
            key = param_name
            is_git_ref = params.get('git-ref', True)
            default = params.get('default', None)
            value_options = params.get('value_options', None)

            parameters.append(Parameter(description, key, value_options, default))

        return cls(parameters)

    def __iter__(self):
        return iter(self.parameters)

    def __getitem__(self, index):
        return self.parameters[index]

    def __setitem__(self, index, value):
        self.parameters[index] = value

    def __len__(self):
        return len(self.parameters)


    def __repr__(self) -> str:
        return f"{self.parameters}"





@dataclass
class Component:
    """
    Represents a component like e.g the openfoam-adapter
    """

    name: str
    template: str
    repository: str
    parameters: Parameters

    def __eq__(self, other):
        if isinstance(other, Component):
            return self.name == other.name
        return False

    def __repr__(self) -> str:
        return f"{self.name}"


class Components(list):
    """
    Represents the collection of components read in from the components.yaml
    """

    def __init__(self, components: List[Component]):
        self.components = components

    @classmethod
    def from_yaml(cls, path):
        """
        Creates a Components instance from a YAML file.

        Args:
            path: The path to the YAML file.

        Returns:
            An instance of Components.
        """
        components = []
        with open(path, 'r') as f:
            data = yaml.safe_load(f)
            for component_name in data:
                parameters = Parameters.from_components_yaml(data[component_name])
                repository = data[component_name]["repository"]
                template = data[component_name]["template"]
                components.append(Component(component_name, template,repository, parameters))

        return cls(components)

    def __iter__(self):
        return iter(self.components)

    def __getitem__(self, index):
        return self.components[index]

    def __setitem__(self, index, value):
        self.components[index] = value

    def __len__(self):
        return len(self.components)

    def get_by_name(self, name_to_search):
        """
        Retrieves a component by its name.

        Args:
            name_to_search: The name of the component to search for.

        Returns:
            The component with the specified name, or None if not found.
        """
        for component in self.components:
            if component.name == name_to_search:
                return component

        return None