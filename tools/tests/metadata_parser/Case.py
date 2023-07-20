from dataclasses import dataclass,field
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from . import Tutorial, Component

@dataclass
class Case:
    """
    Represents a case inside of a tutorial. 
    """

    name: str
    participant: str
    path: Path
    run_cmd: str
    tutorial: Tutorial = field(init=False)
    component: Component

    def __post_init__(self):
        """
        Performs sanity checks after initializing the Case instance.
        """
        if not self.component:
            raise Exception(f'Tried to instantiate the case {self.name} but failed. Reason: Could not find the component it uses in the components.yaml file.')

    @classmethod
    def from_dict(cls, name, dict, available_components):
        """
        Creates a Case instance from a the tutorial yaml dict.

        Args:
            name: The name of the case.
            dict: The dictionary containing the case data.
            available_components: Components read from the components.yaml file

        Returns:
            An instance of the Case but without the tutorial set, this needs to be done later
        """
        participant = dict["participant"]
        path = Path(dict["directory"])
        run_cmd = dict["run"]

        component = available_components.get_by_name(dict["component"])
        return cls(name, participant, path, run_cmd, component)

    def __repr__(self) -> str:
        return f"{self.name}"
