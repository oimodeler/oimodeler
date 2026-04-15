# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 11:01:16 2026

@author: ame
"""

import re
from pathlib import Path

def bump_patch(version):
    major, minor, patch = map(int, version.split("."))
    patch += 1
    return f"{major}.{minor}.{patch}"

def update_file(path, pattern):
    content = Path(path).read_text()

    match = re.search(pattern, content)
    if not match:
        raise ValueError(f"Version non trouvée dans {path}")

    old_version = match.group(1)
    new_version = bump_patch(old_version)

    new_content = re.sub(pattern, f'\\1"{new_version}"', content)
    Path(path).write_text(new_content)

    return old_version, new_version

# fichiers à modifier
files = [
    ("package/__init__.py", r'(__version__\s*=\s*)"([^"]+)"'),
    ("docs/conf.py", r'(release\s*=\s*)"([^"]+)"'),
]

for file, pattern in files:
    old, new = update_file(file, pattern)

print(f"Version mise à jour: {old} -> {new}")