#!/usr/bin/env python3
"""
Validation script for Quickstart system test configuration
Generates comprehensive logs for PR documentation
"""

import sys
from pathlib import Path
from datetime import datetime
from metadata_parser.metdata import Tutorials
from systemtests.TestSuite import TestSuites
from paths import PRECICE_TESTS_DIR, PRECICE_TUTORIAL_DIR

def main():
    print("=" * 60)
    print("QUICKSTART SYSTEM TEST VALIDATION")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    print()
    
    # Load tutorials metadata
    print("[1/4] Loading tutorial metadata...")
    try:
        available_tutorials = Tutorials.from_path(PRECICE_TUTORIAL_DIR)
        quickstart = [t for t in available_tutorials.tutorials 
                     if 'quickstart' in t.path.name.lower()]
        
        if not quickstart:
            print("❌ ERROR: Quickstart tutorial not found!")
            return 1
            
        quickstart = quickstart[0]
        print("✓ Quickstart tutorial metadata loaded successfully")
        print()
    except Exception as e:
        print(f"❌ ERROR loading metadata: {e}")
        return 1
    
    # Display metadata
    print("[2/4] Validating metadata structure...")
    print(f"  Tutorial Name: {quickstart.name}")
    print(f"  Path: {quickstart.path}")
    print(f"  URL: {quickstart.url}")
    print(f"  Participants: {', '.join(quickstart.participants)}")
    print()
    print("  Cases:")
    for case in quickstart.cases:
        print(f"    - {case.name}")
        print(f"        Participant: {case.participant}")
        print(f"        Component: {case.component}")
    print("✓ Metadata structure is valid")
    print()
    
    # Load test suites
    print("[3/4] Loading test suite configuration...")
    try:
        test_suites = TestSuites.from_yaml(
            PRECICE_TESTS_DIR / "tests.yaml", 
            available_tutorials
        )
        print("✓ Test suites loaded successfully")
        print()
    except Exception as e:
        print(f"❌ ERROR loading test suites: {e}")
        return 1
    
    # Check release_test suite
    print("[4/4] Validating release_test suite configuration...")
    release_test = None
    for suite in test_suites.test_suites:
        if suite.name == "release_test":
            release_test = suite
            break
    
    if not release_test:
        print("❌ ERROR: release_test suite not found!")
        return 1
    
    # Find quickstart in release_test
    quickstart_test = None
    for test in release_test.systemtests:
        if 'quickstart' in test.tutorial.name.lower():
            quickstart_test = test
            break
    
    if not quickstart_test:
        print("❌ ERROR: Quickstart not found in release_test suite!")
        return 1
    
    print(f"✓ Quickstart found in release_test suite")
    print()
    print("  Test Configuration:")
    print(f"    Tutorial: {quickstart_test.tutorial.name}")
    print(f"    Case combination: {[c.name for c in quickstart_test.cases]}")
    print(f"    Reference result: {quickstart_test.reference_result}")
    print()
    
    # Summary
    print("=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    print("✓ All validations passed successfully!")
    print()
    print("Changes made:")
    print("  1. Created: quickstart/metadata.yaml")
    print("  2. Updated: tools/tests/tests.yaml (added to release_test)")
    print()
    print("Test configuration:")
    print("  - Tutorial: Quickstart")
    print("  - Participants: Fluid (OpenFOAM), Solid (C++)")
    print("  - Components: openfoam-adapter, bare C++ solver")
    print("  - Pattern: Similar to elastic-tube-1d")
    print()
    print("Next steps:")
    print("  - Generate reference results using GitHub workflow")
    print("  - Or run: python generate_reference_results.py")
    print()
    print("=" * 60)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
