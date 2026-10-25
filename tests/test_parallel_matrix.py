#!/usr/bin/env python3
"""Tests for system test matrix generation."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

from parallel_matrix import (
    build_matrix,
    build_parallel_matrix,
    build_serial_matrix,
    matrix_jobs_for_github,
    total_case_count,
)

TESTS_DIR = Path(__file__).resolve().parent


class ParallelMatrixTests(unittest.TestCase):
    def test_serial_release_is_single_job(self) -> None:
        jobs = build_matrix("release", parallel=False)
        self.assertEqual(len(jobs), 1)
        self.assertEqual(jobs[0].suites, "release")
        self.assertEqual(jobs[0].tutorial, "")

    def test_serial_custom_suites_passthrough(self) -> None:
        jobs = build_matrix("quickstart,elastic-tube-1d", parallel=False)
        self.assertEqual(len(jobs), 1)
        self.assertEqual(jobs[0].suites, "quickstart,elastic-tube-1d")

    def test_parallel_release_splits_by_tutorial(self) -> None:
        jobs = build_matrix("release", parallel=True)
        self.assertGreater(len(jobs), 1)
        tutorials = {job.tutorial for job in jobs}
        self.assertIn("elastic-tube-1d", tutorials)
        self.assertEqual(total_case_count(jobs), 73)

    def test_parallel_quickstart_is_single_tutorial_job(self) -> None:
        jobs = build_matrix("quickstart", parallel=True)
        self.assertEqual(len(jobs), 1)
        self.assertEqual(jobs[0].tutorial, "quickstart")
        self.assertEqual(jobs[0].case_count, 1)

    def test_parallel_elastic_tube_meta_suite(self) -> None:
        jobs = build_matrix("elastic-tube-1d", parallel=True)
        self.assertEqual(len(jobs), 1)
        self.assertEqual(jobs[0].tutorial, "elastic-tube-1d")
        self.assertEqual(jobs[0].case_count, 6)

    def test_parallel_unknown_suite_raises(self) -> None:
        with self.assertRaises(KeyError):
            build_matrix("does-not-exist", parallel=True)

    def test_github_output_is_valid_json(self) -> None:
        jobs = build_parallel_matrix("release")
        payload = matrix_jobs_for_github(matrix_jobs=jobs)
        encoded = json.dumps(payload)
        decoded = json.loads(encoded)
        self.assertEqual(len(decoded), len(jobs))
        self.assertIn("tutorial", decoded[0])
        self.assertIn("suites", decoded[0])

    def test_prepare_matrix_cli_serial(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                "prepare_matrix.py",
                "--suites",
                "release",
                "--mode",
                "serial",
                "--format",
                "github",
            ],
            cwd=TESTS_DIR,
            check=True,
            capture_output=True,
            text=True,
        )
        payload = json.loads(result.stdout.strip())
        self.assertEqual(payload, [{"tutorial": "", "suites": "release"}])

    def test_prepare_matrix_cli_parallel_release(self) -> None:
        result = subprocess.run(
            [
                sys.executable,
                "prepare_matrix.py",
                "--suites",
                "release",
                "--mode",
                "parallel",
                "--format",
                "github",
            ],
            cwd=TESTS_DIR,
            check=True,
            capture_output=True,
            text=True,
        )
        payload = json.loads(result.stdout.strip())
        self.assertGreater(len(payload), 1)
        self.assertEqual(
            sum(len(entry["suites"].split(",")) for entry in payload),
            73,
        )

    def test_build_serial_matrix_helper(self) -> None:
        jobs = build_serial_matrix("extra")
        self.assertEqual(jobs[0].suites, "extra")


if __name__ == "__main__":
    unittest.main()
