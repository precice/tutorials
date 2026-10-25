#!/usr/bin/env python3
"""Emit a GitHub Actions matrix for system test workflows."""

from __future__ import annotations

import argparse
import json
import os
import sys

from parallel_matrix import build_matrix, matrix_jobs_for_github


def _write_github_step_summary(matrix_jobs, suites: str, parallel: bool) -> None:
    summary_path = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary_path:
        return

    job_count = len(matrix_jobs)
    total_cases = sum(job.case_count for job in matrix_jobs)
    with open(summary_path, "w", encoding="utf-8") as summary_file:
        print("## System test matrix\n", file=summary_file)
        if parallel:
            print(
                f"Next, **{job_count}** jobs will run on the "
                f"self-hosted runner(s), split by tutorial.\n",
                file=summary_file,
            )
        else:
            print(
                f"Next, **{job_count}** serial job will run on the "
                f"self-hosted runner(s).\n",
                file=summary_file,
            )
        print(f"- Suites: `{suites}`", file=summary_file)
        print(f"- Mode: `{'parallel' if parallel else 'serial'}`", file=summary_file)
        print(f"- Jobs: **{job_count}**", file=summary_file)
        if parallel:
            print(f"- Total cases: **{total_cases}**\n", file=summary_file)
            print("| # | Tutorial | Cases |", file=summary_file)
            print("| - | --- | ---: |", file=summary_file)
            for index, job in enumerate(matrix_jobs, start=1):
                label = job.tutorial or suites
                print(
                    f"| {index} | `{label}` | {job.case_count} |",
                    file=summary_file,
                )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a system test matrix from tests.yaml.")
    parser.add_argument(
        "--suites",
        "--source-suite",
        dest="suites",
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
        help="Output format (github: JSON for fromJson()).",
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
        print(json.dumps(matrix_jobs_for_github(matrix_jobs=matrix_jobs)))
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
