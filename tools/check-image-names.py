#!python3

import sys
import os
import subprocess


def is_tracked(path: str) -> bool:
    """
    Return True if the given path is tracked by git, False otherwise.

    This guards the hook against local, untracked files that pre-commit
    might pick up when running outside of the normal commit workflow.
    """
    result = subprocess.run(
        ["git", "ls-files", "--error-unmatch", path],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode == 0


problems = False
for file in sys.argv[1:]:
    # Skip files that are not tracked by git
    if not is_tracked(file):
        continue

    parts = file.split(os.sep)
    # Ignore non-interesting files
    if len(parts) != 3 or parts[1] != "images":
        continue

    if parts[0] == "quickstart":
        prefix = "quickstart-"
    else:
        prefix = f"tutorials-{parts[0]}-"

    if not parts[2].startswith(prefix):
        print(f"Incorrect: {file}")
        print(f"Expected prefix: {prefix}")
        print()
        problems = True

sys.exit(1 if problems else 0)
