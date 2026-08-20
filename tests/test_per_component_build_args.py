#!/usr/bin/env python3
"""Unit tests for per-component build arguments (issue #873)."""

from __future__ import annotations

import logging
import os
import unittest
from pathlib import Path

from metadata_parser.metdata import (
    CaseCombination,
    Components,
    ReferenceResult,
    Tutorial,
)
from paths import PRECICE_TESTS_DIR, PRECICE_TUTORIAL_DIR
from systemtests.Systemtest import Systemtest
from systemtests.SystemtestArguments import SystemtestArguments

logging.disable(logging.CRITICAL)


def _patch_param_default(case_combination, component_name: str, key: str, value: str) -> None:
    for case in case_combination.cases:
        if case.component.name != component_name:
            continue
        for param in case.component.parameters:
            if param.key == key:
                param.default = value


def _load_quickstart_systemtest(
    build_args: SystemtestArguments | None = None,
    patches: list[tuple[str, str, str]] | None = None,
) -> tuple[Systemtest, CaseCombination]:
    components = Components.from_yaml(PRECICE_TESTS_DIR / "components.yaml")
    tutorial = Tutorial.from_yaml(
        PRECICE_TUTORIAL_DIR / "quickstart" / "metadata.yaml",
        components,
    )
    case_combination = CaseCombination.from_string_list(
        ["fluid-openfoam", "solid-cpp"], tutorial)
    reference_result = ReferenceResult(
        Path("quickstart/reference-results/fluid-openfoam_solid-cpp.tar.gz"),
        case_combination,
    )

    if patches:
        for component_name, key, value in patches:
            _patch_param_default(case_combination, component_name, key, value)

    systemtest = Systemtest(
        tutorial,
        build_args or SystemtestArguments({}),
        case_combination,
        reference_result,
    )
    systemtest.run_directory = PRECICE_TESTS_DIR.parent / "runs" / "issue873-test"
    systemtest.tutorial_folder = "quickstart_issue873_test"
    return systemtest, case_combination


class PerComponentBuildArgsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        alt_platform = PRECICE_TESTS_DIR / "dockerfiles" / "ubuntu_2404_bare"
        target = PRECICE_TESTS_DIR / "dockerfiles" / "ubuntu_2404"
        if not alt_platform.exists():
            alt_platform.symlink_to(target, target_is_directory=True)

    def test_same_defaults_regression(self) -> None:
        systemtest, _ = _load_quickstart_systemtest()
        bare = systemtest.build_arguments_by_component["bare"]
        openfoam = systemtest.build_arguments_by_component["openfoam-adapter"]
        self.assertEqual(bare.get("PLATFORM"), "ubuntu_2404")
        self.assertEqual(openfoam.get("PLATFORM"), "ubuntu_2404")
        self.assertIsNot(systemtest.params_to_use, systemtest.arguments.arguments)

    def test_different_platform_per_component(self) -> None:
        systemtest, _ = _load_quickstart_systemtest(
            patches=[("bare", "PLATFORM", "ubuntu_2404_bare")])

        bare = systemtest.build_arguments_by_component["bare"]
        openfoam = systemtest.build_arguments_by_component["openfoam-adapter"]
        self.assertEqual(bare.get("PLATFORM"), "ubuntu_2404_bare")
        self.assertEqual(openfoam.get("PLATFORM"), "ubuntu_2404")

        compose = systemtest._Systemtest__get_docker_compose_file()
        self.assertIn("ubuntu_2404_bare", compose)
        self.assertIn("ubuntu_2404", compose.replace("ubuntu_2404_bare", ""))

    def test_different_precice_ref_per_component(self) -> None:
        systemtest, _ = _load_quickstart_systemtest(
            patches=[
                ("bare", "PRECICE_REF", "bare-ref"),
                ("openfoam-adapter", "PRECICE_REF", "openfoam-ref"),
            ])

        bare = systemtest.build_arguments_by_component["bare"]
        openfoam = systemtest.build_arguments_by_component["openfoam-adapter"]
        self.assertEqual(bare.get("PRECICE_REF"), "bare-ref")
        self.assertEqual(openfoam.get("PRECICE_REF"), "openfoam-ref")

        compose = systemtest._Systemtest__get_docker_compose_file()
        self.assertIn("PRECICE_REF=bare-ref", compose)
        self.assertIn("PRECICE_REF=openfoam-ref", compose)

    def test_global_cli_override(self) -> None:
        systemtest, _ = _load_quickstart_systemtest(
            SystemtestArguments({"PLATFORM": "ubuntu_2404"}),
            patches=[("bare", "PLATFORM", "ubuntu_2404_bare")],
        )

        for component_name in ("bare", "openfoam-adapter"):
            platform = systemtest.build_arguments_by_component[component_name].get(
                "PLATFORM")
            self.assertEqual(platform, "ubuntu_2404")


if __name__ == "__main__":
    os.chdir(PRECICE_TESTS_DIR)
    unittest.main()
