#!/usr/bin/env python3
"""
Validate that requirements-reference.txt matches tutorial requirements.txt constraints.

Delegates to report_tutorial_requirements.py --check.
"""
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).parent / "report_tutorial_requirements.py"


def main() -> int:
    result = subprocess.run([sys.executable, str(SCRIPT), "--check"], check=False)
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
