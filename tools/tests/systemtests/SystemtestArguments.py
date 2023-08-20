from dataclasses import dataclass
import yaml
from typing import Optional


@dataclass
class SystemtestArguments:
    arguments: dict[str, str]

    @classmethod
    def from_args(cls, cmd_args):
        if not cmd_args:
            return cls({})

        params_provided = cmd_args.split(",")
        arguments = {}
        for param in params_provided:
            key, value = param.split(":")
            arguments[key] = value

        return cls(arguments)

    @classmethod
    def from_yaml(cls, yml_file):
        if not yml_file:
            return cls({})
        arguments = {}
        with open(yml_file, 'r') as f:
            arguments = yaml.safe_load(f)
        return cls(arguments)

    def __repr__(self):
        return f"{self.arguments}"

    def contains(self, argument_key):
        return argument_key in self.arguments.keys()

    def get(self, argument_key) -> Optional[str]:
        return self.arguments[argument_key]
