#!/usr/bin/env python3
"""
Validate that requirements-reference.txt exists and pyprecice version
matches PYTHON_BINDINGS_REF in reference_versions.yaml.

Exit 0 on success, 1 on failure.
"""
import re
import sys
from pathlib import Path

TOOLS_TESTS = Path(__file__).parent
REFERENCE_VERSIONS = TOOLS_TESTS / "reference_versions.yaml"
REQUIREMENTS_REF = TOOLS_TESTS / "requirements-reference.txt"


def main() -> int:
    if not REQUIREMENTS_REF.exists():
        print(f"ERROR: {REQUIREMENTS_REF} not found. Run update_requirements_reference.py.", file=sys.stderr)
        return 1

    # Load PYTHON_BINDINGS_REF
    ref_text = REFERENCE_VERSIONS.read_text()
    ref_match = re.search(r'PYTHON_BINDINGS_REF:\s*["\']([^"\']+)["\']', ref_text)
    if not ref_match:
        print("ERROR: PYTHON_BINDINGS_REF not found in reference_versions.yaml", file=sys.stderr)
        return 1

    expected_ver = ref_match.group(1).lstrip("v").strip()

    # Parse pyprecice from requirements-reference.txt
    req_text = REQUIREMENTS_REF.read_text()
    precice_match = re.search(r"pyprecice\s*==\s*([\w.]+)", req_text)
    if not precice_match:
        print("ERROR: pyprecice not found in requirements-reference.txt", file=sys.stderr)
        return 1

    actual_ver = precice_match.group(1).strip()
    if actual_ver != expected_ver:
        print(
            f"ERROR: pyprecice version mismatch: requirements-reference.txt has {actual_ver}, "
            f"reference_versions.yaml PYTHON_BINDINGS_REF has {ref_match.group(1)}. "
            "Run: python3 update_requirements_reference.py",
            file=sys.stderr,
        )
        return 1

    print(f"OK: requirements-reference.txt pyprecice=={actual_ver} matches reference_versions.yaml")
    return 0


if __name__ == "__main__":
    sys.exit(main())
