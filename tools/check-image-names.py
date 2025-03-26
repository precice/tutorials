#!python3

import sys
import os

problems = False
for file in sys.argv[1:]:
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
