#!/usr/bin/env python3
"""
Update requirements-reference.txt with pinned versions from reference_versions.yaml.

This script reads PYTHON_BINDINGS_REF from reference_versions.yaml and generates
a requirements-reference.txt with pyprecice pinned to the corresponding version.
Other common packages (numpy, matplotlib, setuptools) use fixed versions known
to work with the tutorials.

Run from tools/tests/:
  python3 update_requirements_reference.py

Or to regenerate from a pip freeze (e.g. from Docker):
  pip freeze | python3 update_requirements_reference.py --from-freeze
"""
import argparse
import re
import sys
from pathlib import Path

REFERENCE_VERSIONS = Path(__file__).parent / "reference_versions.yaml"
REQUIREMENTS_REF = Path(__file__).parent / "requirements-reference.txt"

# Default pinned versions for common packages (fallback when not using --from-freeze)
DEFAULTS = {
    "matplotlib": "3.9.0",
    "numpy": "1.26.4",
    "nutils": "7.2",
    "pyprecice": None,  # From reference_versions.yaml
    "setuptools": ">=69.0.0",
}

# Packages to include (in order)
PACKAGES = ["matplotlib", "numpy", "nutils", "pyprecice", "setuptools"]


def get_pyprecice_version_from_ref(ref: str) -> str:
    """Convert PYTHON_BINDINGS_REF (e.g. v3.2.0) to pyprecice version (3.2.0)."""
    return ref.lstrip("v").strip()


def load_reference_versions() -> str:
    """Load PYTHON_BINDINGS_REF from reference_versions.yaml."""
    text = REFERENCE_VERSIONS.read_text()
    for line in text.splitlines():
        if "PYTHON_BINDINGS_REF" in line and ":" in line:
            match = re.search(r'["\']([^"\']+)["\']', line)
            if match:
                return match.group(1)
    raise ValueError("PYTHON_BINDINGS_REF not found in reference_versions.yaml")


def parse_freezed_packages(freezed: str) -> dict[str, str]:
    """Parse pip freeze output into {package: version}."""
    result = {}
    for line in freezed.strip().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "==" in line:
            pkg, ver = line.split("==", 1)
            result[pkg.lower()] = f"=={ver.strip()}"
        elif "===" in line:
            pkg, ver = line.split("===", 1)
            result[pkg.lower()] = f"=={ver.strip()}"
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Update requirements-reference.txt from reference_versions.yaml"
    )
    parser.add_argument(
        "--from-freeze",
        action="store_true",
        help="Read pip freeze from stdin and use those versions for known packages",
    )
    args = parser.parse_args()

    pyprecice_ref = load_reference_versions()
    pyprecice_ver = get_pyprecice_version_from_ref(pyprecice_ref)

    if args.from_freeze:
        freezed = parse_freezed_packages(sys.stdin.read())
        versions = {}
        for pkg in PACKAGES:
            if pkg.lower() in freezed:
                versions[pkg] = freezed[pkg.lower()]
            elif pkg == "pyprecice":
                versions[pkg] = f"=={pyprecice_ver}"
            elif DEFAULTS.get(pkg):
                versions[pkg] = (
                    DEFAULTS[pkg] if DEFAULTS[pkg].startswith(("==", ">=", "~=")) else f"=={DEFAULTS[pkg]}"
                )
    else:
        DEFAULTS["pyprecice"] = pyprecice_ver
        versions = {
            pkg: f"=={ver}" if ver and not ver.startswith(("==", ">=", "~=")) else (ver or "")
            for pkg, ver in DEFAULTS.items()
        }
        versions["pyprecice"] = f"=={pyprecice_ver}"

    header = """# Pinned Python dependency versions for reproducible system tests and distributions.
# Generated from reference_versions.yaml (PYTHON_BINDINGS_REF). Update at each release.
# Run: python3 update_requirements_reference.py
#
# See tools/tests/README.md for how to update this file.

"""
    lines = [f"{pkg}{versions.get(pkg, '')}\n" for pkg in PACKAGES if pkg in versions]
    REQUIREMENTS_REF.write_text(header + "".join(lines))
    print(f"Wrote {REQUIREMENTS_REF} (pyprecice from {pyprecice_ref})")


if __name__ == "__main__":
    main()
