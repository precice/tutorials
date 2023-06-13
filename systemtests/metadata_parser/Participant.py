from dataclasses import dataclass

@dataclass
class Participant:
    """Represents a participant in a coupled simulation"""

    name: str
    """The name of the participant."""

    def __eq__(self, other) -> bool:
        if isinstance(other, Participant):
            return self.name == other.name
        return False

    def __repr__(self) -> str:
        return f"{self.name}"
