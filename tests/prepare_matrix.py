#!/usr/bin/env python3
"""Build a GitHub Actions matrix for system test workflows.

Reads suite definitions from tests.yaml and emits either:
- serial mode: one job running the requested suite string as-is
- parallel mode: one job per tutorial, with comma-separated per-case suite names

Naming convention for per-case suites:
  <tutorial>_<participant>_<participant>
Example: quickstart_fluid-openfoam_solid-cpp (tutorial: quickstart)
"""

from __future__ import annotations

import argparse
import json
import os
import sys
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


def load_test_suites(tests_yaml: Path | None = None) -> dict:
    path = tests_yaml if tests_yaml is not None else PRECICE_TESTS_DIR / "tests.yaml"
    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    return data["test_suites"]


def build_suite_lookup(test_suites: dict) -> Dict[CaseKey, str]:
    """Map a tutorial case to its runnable suite name in tests.yaml.

    Only considers suite definitions with exactly one tutorial or external
    entry (per-case suites). Meta suites such as release list many cases and
    are skipped here.

    Example:
      ("quickstart", ("fluid-openfoam", "solid-cpp"))
      -> "quickstart_fluid-openfoam_solid-cpp"
    """
    candidates: Dict[CaseKey, set[str]] = defaultdict(set)
    for suite_name, suite_def in test_suites.items():
        for section in ("tutorials", "external"):
            entries = suite_def.get(section, [])
            if len(entries) != 1:
                continue
            entry = entries[0]
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


def _meta_suite_case_entries(
    meta_suite: dict,
    suite_lookup: Dict[CaseKey, str],
) -> List[Tuple[str, str]]:
    """Expand a meta suite (e.g. release) into (tutorial, per-case suite) pairs."""
    entries: List[Tuple[str, str]] = []

    for section in ("tutorials", "external"):
        for case_def in meta_suite.get(section, []):
            case_key: CaseKey = (
                case_def["path"],
                tuple(case_def["case_combination"]),
            )
            tutorial = case_def["path"]
            if case_key not in suite_lookup:
                raise KeyError(
                    f"No runnable suite found for {tutorial} "
                    f"{case_def['case_combination']}")
            entries.append((tutorial, suite_lookup[case_key]))

    return entries


def _tutorial_from_suite_name(suite_name: str) -> str:
    """Derive tutorial name from a per-case suite name.

    Example: channel-transport_fluid-openfoam_transport-openfoam -> channel-transport
    """
    if "_" in suite_name:
        return suite_name.split("_", 1)[0]
    return suite_name


def _is_expandable_meta_suite(test_suites: dict, suite_name: str) -> bool:
    if suite_name not in test_suites:
        return False
    suite_def = test_suites[suite_name]
    return bool(suite_def.get("tutorials") or suite_def.get("external"))


def _build_serial_matrix(suites: str) -> List[MatrixJob]:
    return [MatrixJob(tutorial="", suites=suites, case_count=0)]


def _expand_meta_suite(
    meta_suite_name: str,
    tests_yaml: Path | None = None,
) -> List[MatrixJob]:
    """Split one meta suite into per-tutorial matrix jobs."""
    test_suites = load_test_suites(tests_yaml)
    if meta_suite_name not in test_suites:
        raise KeyError(f"Unknown test suite: {meta_suite_name}")

    suite_lookup = build_suite_lookup(test_suites)
    case_entries = _meta_suite_case_entries(
        test_suites[meta_suite_name],
        suite_lookup,
    )

    grouped: Dict[str, List[str]] = defaultdict(list)
    for tutorial, case_suite in case_entries:
        grouped[tutorial].append(case_suite)

    return [
        MatrixJob(
            tutorial=tutorial,
            suites=",".join(case_suites),
            case_count=len(case_suites),
        )
        for tutorial in sorted(grouped)
        for case_suites in [grouped[tutorial]]
    ]


def _merge_matrix_jobs(jobs: List[MatrixJob]) -> List[MatrixJob]:
    """Merge jobs for the same tutorial when multiple meta suites were requested."""
    grouped_suites: Dict[str, List[str]] = defaultdict(list)
    for job in jobs:
        for suite_name in job.suites.split(","):
            if suite_name:
                grouped_suites[job.tutorial].append(suite_name)

    merged: List[MatrixJob] = []
    for tutorial in sorted(grouped_suites):
        unique_suites: List[str] = []
        seen: set[str] = set()
        for suite_name in grouped_suites[tutorial]:
            if suite_name in seen:
                continue
            seen.add(suite_name)
            unique_suites.append(suite_name)
        merged.append(
            MatrixJob(
                tutorial=tutorial,
                suites=",".join(unique_suites),
                case_count=len(unique_suites),
            )
        )
    return merged


def build_matrix(
    suites: str,
    parallel: bool,
    tests_yaml: Path | None = None,
) -> List[MatrixJob]:
    if not parallel:
        return _build_serial_matrix(suites)

    test_suites = load_test_suites(tests_yaml)
    suite_names = [part.strip() for part in suites.split(",") if part.strip()]
    if not suite_names:
        raise KeyError("No test suites specified")

    matrix_jobs: List[MatrixJob] = []
    for suite_name in suite_names:
        if _is_expandable_meta_suite(test_suites, suite_name):
            matrix_jobs.extend(_expand_meta_suite(suite_name, tests_yaml))
            continue
        if suite_name in test_suites:
            matrix_jobs.append(
                MatrixJob(
                    tutorial=_tutorial_from_suite_name(suite_name),
                    suites=suite_name,
                    case_count=1,
                )
            )
            continue
        raise KeyError(f"Unknown test suite: {suite_name}")

    if len(suite_names) == 1:
        return matrix_jobs
    return _merge_matrix_jobs(matrix_jobs)


def matrix_jobs_for_github(matrix_jobs: List[MatrixJob]) -> List[dict]:
    return [{"tutorial": job.tutorial, "suites": job.suites} for job in matrix_jobs]


def _write_github_step_summary(
    matrix_jobs: List[MatrixJob],
    suites: str,
    parallel: bool,
) -> None:
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary_path:
        return

    job_count = len(matrix_jobs)
    with open(summary_path, "w", encoding="utf-8") as summary_file:
        print("## System test matrix\n", file=summary_file)
        if parallel:
            print(
                f"Next, **{job_count}** jobs will run, split by tutorial.\n",
                file=summary_file,
            )
        else:
            assert job_count == 1
            print(
                "Next, all test cases will be executed as one serial job.\n",
                file=summary_file,
            )
        print(f"- Suites: `{suites}`", file=summary_file)
        print(f"- Mode: `{'parallel' if parallel else 'serial'}`", file=summary_file)
        print(f"- Jobs: **{job_count}**", file=summary_file)
        if parallel:
            total_cases = sum(job.case_count for job in matrix_jobs)
            print(f"- Total cases: **{total_cases}**\n", file=summary_file)
            print("| # | Tutorial | Cases |", file=summary_file)
            print("| - | --- | ---: |", file=summary_file)
            for index, job in enumerate(matrix_jobs, start=1):
                print(
                    f"| {index} | `{job.tutorial}` | {job.case_count} |",
                    file=summary_file,
                )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a system test matrix from tests.yaml.")
    parser.add_argument(
        "-s",
        "--suites",
        default="release",
        help="Suite or comma-separated suites to run (default: release).",
    )
    parser.add_argument(
        "--mode",
        choices=["serial", "parallel"],
        default="parallel",
        help="serial: one matrix job; parallel: split expandable suites per tutorial.",
    )
    parser.add_argument(
        "--format",
        choices=["github", "pretty"],
        default="pretty",
        help=(
            "github: JSON for workflow fromJson(); "
            "pretty: human-readable output for local runs."
        ),
    )
    args = parser.parse_args()

    parallel = args.mode == "parallel"
    try:
        matrix_jobs = build_matrix(args.suites, parallel=parallel)
    except (KeyError, FileNotFoundError) as error:
        print(error, file=sys.stderr)
        return 1

    _write_github_step_summary(matrix_jobs, args.suites, parallel)

    if args.format == "github":
        print(json.dumps(matrix_jobs_for_github(matrix_jobs)))
    else:
        total_cases = sum(job.case_count for job in matrix_jobs)
        print(
            f"Matrix for {args.suites} ({args.mode}): "
            f"{len(matrix_jobs)} job(s), {total_cases} case(s)")
        for index, job in enumerate(matrix_jobs, start=1):
            label = job.tutorial or "(serial)"
            print(f"  {index}. {label}: {job.case_count} case(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
