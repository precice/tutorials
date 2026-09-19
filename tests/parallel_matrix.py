"""Build a per-tutorial parallel matrix from a meta test suite in tests.yaml."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import yaml

from paths import PRECICE_TESTS_DIR

CaseKey = Tuple[str, Tuple[str, ...]]


@dataclass(frozen=True)
class MatrixJob:
    tutorial: str
    suites: str
    case_count: int


def _tests_yaml_path(tests_yaml: Path | None = None) -> Path:
    return tests_yaml if tests_yaml is not None else PRECICE_TESTS_DIR / "tests.yaml"


def load_test_suites(tests_yaml: Path | None = None) -> dict:
    with open(_tests_yaml_path(tests_yaml), "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    return data["test_suites"]


def build_case_suite_lookup(test_suites: dict) -> Dict[CaseKey, str]:
    """Map (path, case_combination) to the runnable per-case suite name."""
    candidates: Dict[CaseKey, set[str]] = defaultdict(set)
    for suite_name, suite_def in test_suites.items():
        tutorials = suite_def.get("tutorials", [])
        if len(tutorials) != 1:
            continue
        entry = tutorials[0]
        key = (entry["path"], tuple(entry["case_combination"]))
        candidates[key].add(suite_name)

    lookup: Dict[CaseKey, str] = {}
    for key, suite_names in candidates.items():
        path, combo = key
        canonical = f"{path}_{'_'.join(combo)}"
        if canonical in suite_names:
            lookup[key] = canonical
        else:
            lookup[key] = sorted(suite_names)[0]
    return lookup


def build_external_suite_lookup(test_suites: dict) -> Dict[CaseKey, str]:
    """Map external cases to the suite that defines them (e.g. mixedbc)."""
    lookup: Dict[CaseKey, str] = {}
    for suite_name, suite_def in test_suites.items():
        external = suite_def.get("external", [])
        if len(external) != 1:
            continue
        entry = external[0]
        lookup[(entry["path"], tuple(entry["case_combination"]))] = suite_name
    return lookup


def _release_case_entries(
    source_suite: dict,
    case_lookup: Dict[CaseKey, str],
    external_lookup: Dict[CaseKey, str],
) -> List[Tuple[str, str]]:
    """Return (tutorial_path, case_suite_name) pairs for a meta suite definition."""
    entries: List[Tuple[str, str]] = []

    def add_entry(case_def: dict) -> None:
        case_key: CaseKey = (case_def["path"], tuple(case_def["case_combination"]))
        tutorial = case_def["path"]
        if case_key in case_lookup:
            suite_name = case_lookup[case_key]
        elif case_key in external_lookup:
            suite_name = external_lookup[case_key]
        else:
            raise KeyError(
                f"No runnable suite found for {tutorial} "
                f"{case_def['case_combination']}")
        entries.append((tutorial, suite_name))

    for case_def in source_suite.get("tutorials", []):
        add_entry(case_def)
    for case_def in source_suite.get("external", []):
        add_entry(case_def)

    return entries


def build_parallel_matrix(
    source_suite_name: str = "release",
    tests_yaml: Path | None = None,
) -> List[MatrixJob]:
    test_suites = load_test_suites(tests_yaml)
    if source_suite_name not in test_suites:
        raise KeyError(f"Unknown test suite: {source_suite_name}")

    case_lookup = build_case_suite_lookup(test_suites)
    external_lookup = build_external_suite_lookup(test_suites)
    release_entries = _release_case_entries(
        test_suites[source_suite_name],
        case_lookup,
        external_lookup,
    )

    grouped: Dict[str, List[str]] = defaultdict(list)
    for tutorial, case_suite in release_entries:
        grouped[tutorial].append(case_suite)

    matrix_jobs: List[MatrixJob] = []
    for tutorial in sorted(grouped):
        case_suites = grouped[tutorial]
        matrix_jobs.append(
            MatrixJob(
                tutorial=tutorial,
                suites=",".join(case_suites),
                case_count=len(case_suites),
            )
        )
    return matrix_jobs


def matrix_jobs_for_github(
    source_suite_name: str = "release",
    tests_yaml: Path | None = None,
    matrix_jobs: List[MatrixJob] | None = None,
) -> List[dict]:
    jobs = (
        matrix_jobs
        if matrix_jobs is not None
        else build_parallel_matrix(source_suite_name, tests_yaml))
    return [
        {"tutorial": job.tutorial, "suites": job.suites}
        for job in jobs
    ]
