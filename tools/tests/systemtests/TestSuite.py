from dataclasses import dataclass, field
from typing import Optional, List, Dict
from metadata_parser.metdata import Tutorials, Tutorial, Case, CaseCombination, ReferenceResult

import yaml


@dataclass
class TestSuite:
    name: str
    cases_of_tutorial: Dict[Tutorial, List[CaseCombination]]
    reference_results: Dict[Tutorial, List[ReferenceResult]]

    def __repr__(self) -> str:
        return_string = f"Test suite: {self.name} contains:"
        for tutorial, cases in self.cases_of_tutorial.items():
            return_string += f"""
    {tutorial.name}
        cases: {cases}
        reference_results: {self.reference_results[tutorial]}"""

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
            test_suites_raw = data['test_suites']
            for test_suite_name in test_suites_raw:
                case_combinations_of_tutorial = {}
                reference_results_of_tutorial = {}
                # iterate over tutorials:
                for tutorial_case in test_suites_raw[test_suite_name]['tutorials']:
                    tutorial = parsed_tutorials.get_by_path(tutorial_case['path'])
                    if not tutorial:
                        raise Exception(f"No tutorial with path {tutorial_case['path']} found.")
                    # initialize the datastructure for the new Testsuite
                    if tutorial not in case_combinations_of_tutorial:
                        case_combinations_of_tutorial[tutorial] = []
                        reference_results_of_tutorial[tutorial] = []

                    all_case_combinations = tutorial.case_combinations
                    case_combination_requested = CaseCombination.from_string_list(
                        tutorial_case['case_combination'], tutorial)
                    if case_combination_requested in all_case_combinations:
                        case_combinations_of_tutorial[tutorial].append(case_combination_requested)
                        reference_results_of_tutorial[tutorial].append(ReferenceResult(
                            tutorial_case['reference_result'], case_combination_requested))
                    else:
                        raise Exception(
                            f"Could not find the following cases {tutorial_case['case-combination']} in the current metadata of tutorial {tutorial.name}")

                testsuites.append(TestSuite(test_suite_name, case_combinations_of_tutorial,
                                            reference_results_of_tutorial))

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
