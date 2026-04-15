#!/usr/bin/env python3

import re
from pathlib import Path


def bump_patch(v):
    major, minor, patch = map(int, v.split("."))
    return f"{major}.{minor}.{patch + 1}"


path = Path("oimodeler/__init__.py")
content = path.read_text()

match = re.search(r'__version__\s*=\s*"([^"]+)"', content)
if not match:
    raise ValueError("Version not found")

old = match.group(1)
new = bump_patch(old)

content = re.sub(
    r'(__version__\s*=\s*)"[^"]+"',
    rf'\1"{new}"',
    content
)

path.write_text(content)

print(f"{old} → {new}")