import subprocess
import re
from pathlib import Path


PACKAGE_FILE = Path("oimodeler/__init__.py")


def get_commits_since_last_tag():
    """Récupère les messages de commit depuis le dernier tag."""
    try:
        last_tag = subprocess.check_output(
            ["git", "describe", "--tags", "--abbrev=0"],
            text=True
        ).strip()
        cmd = ["git", "log", f"{last_tag}..HEAD", "--pretty=%B"]
    except subprocess.CalledProcessError:
        cmd = ["git", "log", "--pretty=%B"]

    output = subprocess.check_output(cmd, text=True)
    return output.splitlines()


def detect_bump_level(commits):
    """
    Règles :
    - MAJOR: → MAJOR
    - MINOR: → MINOR
    - PATCH: → PATCH
    - sinon → None (pas de bump)
    """
    level = None

    for msg in commits:
        msg = msg.strip()

        if msg.startswith("MAJOR:"):
            return "MAJOR"
        elif msg.startswith("MINOR:"):
            level = "MINOR"
        elif msg.startswith("PATCH:") and level is None:
            level = "PATCH"

    return level


def bump_version(version, level):
    major, minor, patch = map(int, version.split("."))

    if level == "MAJOR":
        return f"{major + 1}.0.0"
    elif level == "MINOR":
        return f"{major}.{minor + 1}.0"
    elif level == "PATCH":
        return f"{major}.{minor}.{patch + 1}"
    else:
        return None


def update_file(old_version, new_version):
    content = PACKAGE_FILE.read_text()

    content = re.sub(
        r'(__version__\s*=\s*)"[^"]+"',
        rf'\1"{new_version}"',
        content
    )

    PACKAGE_FILE.write_text(content)


def main():
    content = PACKAGE_FILE.read_text()

    match = re.search(r'__version__\s*=\s*"([^"]+)"', content)
    if not match:
        raise ValueError("Version not found in __init__.py")

    current_version = match.group(1)

    commits = get_commits_since_last_tag()
    level = detect_bump_level(commits)

    if level is None:
        print("No bump keyword found → skipping version bump")
        return

    new_version = bump_version(current_version, level)

    if new_version is None:
        print("No valid bump → skipping")
        return

    update_file(current_version, new_version)

    print(f"{current_version} → {new_version} ({level})")


if __name__ == "__main__":
    main()