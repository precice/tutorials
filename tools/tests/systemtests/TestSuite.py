from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict
from metadata_parser.metdata import (
    Tutorials,
    Tutorial,
    Case,
    CaseCombination,
    ReferenceResult,
    Components,
)
from paths import PRECICE_TESTS_DIR
from systemtests.sources import (
    TutorialSource,
    resolve_tutorial_root,
    PRECICE_EXTERNAL_CACHE_DIR,
)

import yaml


@dataclass
class TestSuite:
    name: str
    cases_of_tutorial: Dict[Tutorial, List[CaseCombination]]
    reference_results: Dict[Tutorial, List[ReferenceResult]]
    max_times: Dict[Tutorial, list] = field(default_factory=dict)
    max_time_windows: Dict[Tutorial, list] = field(default_factory=dict)
    timeouts: Dict[Tutorial, List] = field(default_factory=dict)
    tolerances: Dict[Tutorial, list] = field(default_factory=dict)
    skip_compares: Dict[Tutorial, list] = field(default_factory=dict)
    run_befores: Dict[Tutorial, List] = field(default_factory=dict)
    run_afters: Dict[Tutorial, List] = field(default_factory=dict)

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
        available_components = Components.from_yaml(PRECICE_TESTS_DIR / "components.yaml")
        with open(path, 'r') as f:
            data = yaml.safe_load(f)
            test_suites_raw = data['test_suites']
            for test_suite_name in test_suites_raw:
                case_combinations_of_tutorial = {}
                reference_results_of_tutorial = {}
                max_times_of_tutorial = {}
                max_time_windows_of_tutorial = {}
                timeouts_of_tutorial = {}
                tolerances_of_tutorial = {}
                skip_compares_of_tutorial = {}
                run_befores_of_tutorial = {}
                run_afters_of_tutorial = {}
                suite_def = test_suites_raw[test_suite_name]
                local_cases = suite_def.get('tutorials', [])
                external_cases = suite_def.get('external', [])
                if not local_cases and not external_cases:
                    raise ValueError(
                        f"Test suite '{test_suite_name}' must define 'tutorials' or 'external'."
                    )

                for tutorial_case in local_cases:
                    cls._add_tutorial_case(
                        tutorial_case,
                        TutorialSource.local(),
                        parsed_tutorials,
                        available_components,
                        case_combinations_of_tutorial,
                        reference_results_of_tutorial,
                        max_times_of_tutorial,
                        max_time_windows_of_tutorial,
                        timeouts_of_tutorial,
                        tolerances_of_tutorial,
                        skip_compares_of_tutorial,
                        run_befores_of_tutorial,
                        run_afters_of_tutorial,
                    )

                for tutorial_case in external_cases:
                    source_raw = tutorial_case.get('source')
                    if not source_raw:
                        raise ValueError(
                            f"External test entry in suite '{test_suite_name}' "
                            f"requires a 'source' block."
                        )
                    source = TutorialSource.from_dict(source_raw)
                    if source.type == "local":
                        raise ValueError(
                            f"External test entry in suite '{test_suite_name}' "
                            f"must use a git or archive source."
                        )
                    cls._add_tutorial_case(
                        tutorial_case,
                        source,
                        parsed_tutorials,
                        available_components,
                        case_combinations_of_tutorial,
                        reference_results_of_tutorial,
                        max_times_of_tutorial,
                        max_time_windows_of_tutorial,
                        timeouts_of_tutorial,
                        tolerances_of_tutorial,
                        skip_compares_of_tutorial,
                        run_befores_of_tutorial,
                        run_afters_of_tutorial,
                    )

                testsuites.append(TestSuite(
                    test_suite_name,
                    case_combinations_of_tutorial,
                    reference_results_of_tutorial,
                    max_times_of_tutorial,
                    max_time_windows_of_tutorial,
                    timeouts_of_tutorial,
                    tolerances_of_tutorial,
                    skip_compares_of_tutorial,
                    run_befores_of_tutorial,
                    run_afters_of_tutorial,
                ))

        return cls(testsuites)

    @staticmethod
    def _add_tutorial_case(
        tutorial_case: dict,
        source: TutorialSource,
        parsed_tutorials: Tutorials,
        available_components: Components,
        case_combinations_of_tutorial: Dict[Tutorial, List[CaseCombination]],
        reference_results_of_tutorial: Dict[Tutorial, List[ReferenceResult]],
        max_times_of_tutorial: Dict[Tutorial, list],
        max_time_windows_of_tutorial: Dict[Tutorial, list],
        timeouts_of_tutorial: Dict[Tutorial, List],
        tolerances_of_tutorial: Dict[Tutorial, list],
        skip_compares_of_tutorial: Dict[Tutorial, list],
        run_befores_of_tutorial: Dict[Tutorial, list],
        run_afters_of_tutorial: Dict[Tutorial, list],
    ) -> None:
        tutorial = parsed_tutorials.get_by_path(tutorial_case['path'])
        if not tutorial and source.type != "local":
            tutorial_root = resolve_tutorial_root(
                Path(tutorial_case['path']),
                source,
                PRECICE_EXTERNAL_CACHE_DIR,
            )
            metadata_path = tutorial_root / "metadata.yaml"
            if not metadata_path.exists():
                raise FileNotFoundError(
                    f"No metadata.yaml found for external tutorial "
                    f"{tutorial_case['path']} at {tutorial_root}"
                )
            tutorial = Tutorial.from_yaml(
                metadata_path,
                available_components,
                base_dir=tutorial_root.parent,
                source=source,
            )
            tutorial.resolved_root = tutorial_root
            parsed_tutorials.tutorials.append(tutorial)
        if not tutorial:
            raise Exception(f"No tutorial with path {tutorial_case['path']} found.")
        if tutorial not in case_combinations_of_tutorial:
            case_combinations_of_tutorial[tutorial] = []
            reference_results_of_tutorial[tutorial] = []
            max_times_of_tutorial[tutorial] = []
            max_time_windows_of_tutorial[tutorial] = []
            timeouts_of_tutorial[tutorial] = []
            tolerances_of_tutorial[tutorial] = []
            skip_compares_of_tutorial[tutorial] = []
            run_befores_of_tutorial[tutorial] = []
            run_afters_of_tutorial[tutorial] = []

        all_case_combinations = tutorial.case_combinations
        case_combination_requested = CaseCombination.from_string_list(
            tutorial_case['case_combination'], tutorial)
        if case_combination_requested in all_case_combinations:
            case_combinations_of_tutorial[tutorial].append(case_combination_requested)
            ref_base = tutorial.path.parent if source.type != "local" else None
            reference_results_of_tutorial[tutorial].append(ReferenceResult(
                tutorial_case['reference_result'],
                case_combination_requested,
                base_dir=ref_base,
            ))
            max_time_raw = tutorial_case.get('max_time', None)
            if max_time_raw is not None and (not isinstance(
                    max_time_raw, (int, float)) or max_time_raw <= 0):
                raise ValueError(f"max_time must be a positive number, got {max_time_raw!r}")
            max_times_of_tutorial[tutorial].append(max_time_raw)
            mtw_raw = tutorial_case.get('max_time_windows', None)
            if mtw_raw is not None and (not isinstance(mtw_raw, int) or mtw_raw <= 0):
                raise ValueError(f"max_time_windows must be a positive integer, got {mtw_raw!r}")
            max_time_windows_of_tutorial[tutorial].append(mtw_raw)

            timeout_value = tutorial_case.get('timeout', None)
            if timeout_value is not None and not isinstance(timeout_value, int):
                raise TypeError(
                    f"Expected 'timeout' to be an integer or None, but got {type(timeout_value).__name__} "
                    f"(value: {timeout_value}) in tutorial '{tutorial}'."
                )
            timeouts_of_tutorial[tutorial].append(timeout_value)

            tolerance_value = tutorial_case.get('tolerance', None)
            if tolerance_value is not None:
                if isinstance(tolerance_value, str):
                    try:
                        tolerance_value = float(tolerance_value)
                    except ValueError as exc:
                        raise ValueError(
                            f"tolerance must be a positive number, got {tolerance_value!r}") from exc
                if not isinstance(tolerance_value, (int, float)) or tolerance_value <= 0:
                    raise ValueError(
                        f"tolerance must be a positive number, got {tolerance_value!r}")
            tolerances_of_tutorial[tutorial].append(tolerance_value)

            skip_compare_value = tutorial_case.get('skip_compare', None)
            if skip_compare_value is not None and not isinstance(skip_compare_value, bool):
                raise TypeError(
                    f"Expected 'skip_compare' to be a boolean or None, but got "
                    f"{type(skip_compare_value).__name__} (value: {skip_compare_value}) "
                    f"in tutorial '{tutorial}'."
                )
            skip_compares_of_tutorial[tutorial].append(skip_compare_value)

            run_before_raw = tutorial_case.get('run-before', None)
            run_after_raw = tutorial_case.get('run-after', None)
            run_befores_of_tutorial[tutorial].append(
                run_before_raw.strip()
                if isinstance(run_before_raw, str) and run_before_raw.strip() else None)
            run_afters_of_tutorial[tutorial].append(
                run_after_raw.strip()
                if isinstance(run_after_raw, str) and run_after_raw.strip() else None)
        else:
            raise Exception(
                f"Could not find the case combination {tutorial_case['case_combination']} "
                f"in the current metadata of tutorial {tutorial.name}, or it does not "
                f"define all necessary participants.")

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
