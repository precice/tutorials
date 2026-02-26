#!/usr/bin/env python3
"""Check components configuration for quickstart"""
import yaml
from paths import PRECICE_TESTS_DIR

print("=" * 60)
print("COMPONENT CONFIGURATION CHECK")
print("=" * 60)

with open(PRECICE_TESTS_DIR / "components.yaml", 'r') as f:
    components = yaml.safe_load(f)

print("\nComponents Required for Quickstart:")
print("-" * 60)

# Quickstart uses openfoam-adapter and bare C++ compiler
required = ["openfoam-adapter", "bare"]

for comp_name in required:
    if comp_name in components:
        print(f"\n✓ {comp_name}")
        comp = components[comp_name]
        if 'repository' in comp:
            print(f"  Repository: {comp['repository']}")
        if 'template' in comp:
            print(f"  Template: {comp['template']}")
    else:
        print(f"\n✗ {comp_name} NOT FOUND")

print("\n" + "=" * 60)
print("Component check complete")
print("=" * 60)
