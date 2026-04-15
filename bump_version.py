#!/usr/bin/env python3

import re
from pathlib import Path


def bump_patch(version: str) -> str:
    major, minor, patch = map(int, version.split("."))
    patch += 1
    return f"{major}.{minor}.{patch}"


def update_file(path: Path, pattern: str):
    if not path.exists():
        raise FileNotFoundError(f"Fichier introuvable: {path}")

    content = path.read_text()

    match = re.search(pattern, content)
    if not match:
        raise ValueError(f"Version non trouvée dans {path}")

    old_version = match.group(2)
    new_version = bump_patch(old_version)

    new_content = re.sub(pattern, rf'\1"{new_version}"', content)
    path.write_text(new_content)

    print(f"{path} : {old_version} → {new_version}")

    return old_version, new_version


def main():
    BASE_DIR = Path(__file__).resolve().parent

    files = [
        (BASE_DIR / "oimodeler" / "__init__.py",
         r'(__version__\s*=\s*)"([^"]+)"'),
    ]

    for path, pattern in files:
        update_file(path, pattern)


if __name__ == "__main__":
    print("=== Bump version PATCH ===")
    main()