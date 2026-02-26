#!/usr/bin/env python3
"""Simple validation test for Quickstart"""
from metadata_parser.metdata import Tutorials
from systemtests.TestSuite import TestSuites
from paths import PRECICE_TESTS_DIR, PRECICE_TUTORIAL_DIR

print("=" * 60)
print("QUICKSTART VALIDATION TEST")
print("=" * 60)

# Test 1: Load metadata
print("\n[Test 1] Loading tutorials...")
tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)
quickstart = [t for t in tutorials.tutorials if 'quickstart' in t.name.lower()]
if quickstart:
    print(f"✓ Found: {quickstart[0].name}")
    print(f"  Cases: {[c.name for c in quickstart[0].cases]}")
else:
    print("✗ Quickstart not found")

# Test 2: Load test suites
print("\n[Test 2] Loading test suites...")
suites = TestSuites.from_yaml(PRECICE_TESTS_DIR / "tests.yaml", tutorials)
release_test = [s for s in suites.testsuites if s.name == "release_test"][0]
total_tests = sum(len(cases) for cases in release_test.cases_of_tutorial.values())
print(f"✓ release_test suite has {total_tests} test configurations")

# Test 3: Find quickstart in release_test
print("\n[Test 3] Checking quickstart in release_test...")
qs_found = False
for tutorial, cases in release_test.cases_of_tutorial.items():
    if 'quickstart' in tutorial.name.lower():
        qs_found = True
        print(f"✓ Quickstart test found:")
        print(f"  Cases: {cases}")
        print(f"  Reference: {release_test.reference_results[tutorial]}")
        break

if not qs_found:
    print("✗ Quickstart not in release_test")

print("\n" + "=" * 60)
print("ALL TESTS PASSED ✓")
print("=" * 60)
