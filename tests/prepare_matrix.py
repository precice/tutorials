#!/usr/bin/env python3
"""Emit a GitHub Actions matrix for parallel release system tests."""

from __future__ import annotations

import argparse
import json
import os
import sys

from parallel_matrix import build_parallel_matrix, matrix_jobs_for_github


def _write_github_step_summary(matrix_jobs) -> None:
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary_path:
        return

    total_cases = sum(job.case_count for job in matrix_jobs)
    with open(summary_path, "w", encoding="utf-8") as summary_file:
        print("## Parallel test matrix\n", file=summary_file)
        print(f"- Source suite: `release`", file=summary_file)
        print(f"- Tutorial jobs: **{len(matrix_jobs)}**", file=summary_file)
        print(f"- Total cases: **{total_cases}**\n", file=summary_file)
        print("| Tutorial | Cases |", file=summary_file)
        print("| --- | ---: |", file=summary_file)
        for job in matrix_jobs:
            print(f"| `{job.tutorial}` | {job.case_count} |", file=summary_file)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a per-tutorial matrix from a meta test suite.")
    parser.add_argument(
        "--source-suite",
        default="release",
        help="Meta suite to expand (default: release).",
    )
    parser.add_argument(
        "--format",
        choices=["github", "pretty"],
        default="pretty",
        help="Output format (github: JSON for fromJson()).",
    )
    args = parser.parse_args()

    try:
        matrix_jobs = build_parallel_matrix(args.source_suite)
    except (KeyError, FileNotFoundError) as error:
        print(error, file=sys.stderr)
        return 1

    _write_github_step_summary(matrix_jobs)

    if args.format == "github":
        print(json.dumps(matrix_jobs_for_github(args.source_suite)))
    else:
        total_cases = sum(job.case_count for job in matrix_jobs)
        print(
            f"Matrix for {args.source_suite}: "
            f"{len(matrix_jobs)} tutorial jobs, {total_cases} cases")
        for job in matrix_jobs:
            print(f"  {job.tutorial}: {job.case_count} case(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
