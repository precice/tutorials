from dataclasses import dataclass
from typing import Optional,List
import yaml

@dataclass
class BuildArgument:
    """Represents a BuildArgument needed to run the docker container"""

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
        Check if the BuildArgument need to be supplied via CommandLineArgs

        Returns:
            bool: True if the parameter is required, False otherwise.
        """
        return False if self.default else True
    
    def __eq__(self, other) -> bool:
        if isinstance(other, BuildArgument):
            return self.key == other.key
        return False

    def __hash__(self) -> int:
        return hash(self.key)

    def __repr__(self) -> str:
        return f"{self.key}"



class BuildArguments:
    """Represents a collection of parameters."""

    def __init__(self, arguments: List[BuildArgument]):
        self.arguments  = arguments

    @classmethod
    def from_components_yaml(cls, data):
        """
        Create a list of Paramters from the components YAML data.

        Args:
            data: The components YAML data.
        """
        arguments = []
        for param_name, params in data['build-arguments'].items():
            ## TODO maybe **params
            description = params.get('description', f"No description provided for {param_name}")
            key = param_name
            default = params.get('default', None)
            value_options = params.get('value_options', None)

            arguments.append(BuildArgument(description, key, value_options, default))

        return cls(arguments)

    def __iter__(self):
        return iter(self.arguments)

    def __getitem__(self, index):
        return self.arguments[index]

    def __setitem__(self, index, value):
        self.arguments[index] = value

    def __len__(self):
        return len(self.arguments)


    def __repr__(self) -> str:
        return f"{self.arguments}"





@dataclass
class Component:
    """
    Represents a component like e.g the openfoam-adapter
    """

    name: str
    template: str
    repository: str
    parameters: BuildArguments

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
                parameters = BuildArguments.from_components_yaml(data[component_name])
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