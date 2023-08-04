from dataclasses import dataclass


@dataclass
class CmdLineArguments:
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

    def __repr__(self):
        return f"{self.arguments}"

    def contains(self, argument_key):
        return argument_key in self.arguments.keys()

    def get(self, argument_key):
        return self.arguments[argument_key]
