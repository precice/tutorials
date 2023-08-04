from dataclasses import dataclass, field
from typing import Optional, List, Dict
from metadata_parser.metdata import Tutorials, Tutorial, Case

import yaml


@dataclass
class TestSuite:
    name: str
    cases_of_tutorial: Dict[Tutorial, List[Case]]

    def __repr__(self) -> str:
        return_string = f"Test suite: {self.name} contains:"
        for tutorial, cases in self.cases_of_tutorial.items():
            return_string += f"""
    {tutorial.name}
        cases: {cases}"""

        return return_string


class TestSuites(list):
    """
    Represents the collection of testsuites read in from the tests.yaml
    """

    def __init__(self, testsuites: List[TestSuite]):
        self.testsuites = testsuites

    @classmethod
    def from_yaml(cls, path, parsed_tutorials: Tutorials):
        """
        Creates a TestSuites instance from a YAML file.

        Args:
            path: The path to the YAML file.

        Returns:
            An instance of TestSuites.
        """
        testsuites = []
        with open(path, 'r') as f:
            data = yaml.safe_load(f)
            test_suites_raw = data['test-suites']
            for test_suite_name in test_suites_raw:
                cases_of_tutorial = {}
                # iterate over tutorials:
                for tutorial_path in test_suites_raw[test_suite_name]['tutorials']:
                    tutorial = parsed_tutorials.get_by_path(tutorial_path)
                    if not tutorial:
                        raise Exception(
                            f"No tutorial with path {tutorial_path} found.")
                    cases_of_tutorial[tutorial] = []
                    all_cases = tutorial.case_combinations
                    cases_requested = test_suites_raw[test_suite_name]['tutorials'][tutorial_path]['cases']
                    for case_combination in all_cases:
                        if f"{case_combination}" in cases_requested:
                            cases_of_tutorial[tutorial].append(
                                case_combination)
                testsuites.append(
                    TestSuite(test_suite_name, cases_of_tutorial))

        return cls(testsuites)

    def __iter__(self):
        return iter(self.testsuites)

    def __getitem__(self, index):
        return self.testsuites[index]

    def __setitem__(self, index, value):
        self.testsuites[index] = value

    def __len__(self):
        return len(self.testsuites)

    def get_by_name(self, name_to_search) -> Optional[TestSuite]:
        """
        Retrieves a testsuite by its name.

        Args:
            name_to_search: The name of the testsuite to search for.

        Returns:
            The component with the specified name, or None if not found.
        """
        for testsuite in self.testsuites:
            if testsuite.name == name_to_search:
                return testsuite

        return None

    def __repr__(self) -> str:
        return_str = ""
        for tests_suite in self.testsuites:
            return_str += f"{tests_suite}\n\n"
        return return_str
