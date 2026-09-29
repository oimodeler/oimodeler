import os
import re
import subprocess
from pathlib import Path


BASE_DIR = Path(__file__).resolve().parent
PACKAGE_FILE = BASE_DIR / "oimodeler" / "__init__.py"


def get_commits():
    """
    Get commits included in the current GitHub push.

    Falls back to commits since the latest tag when running locally.
    """

    github_before = os.environ.get("GITHUB_EVENT_BEFORE")
    github_sha = os.environ.get("GITHUB_SHA")

    # Running in GitHub Actions
    if github_before and github_sha:
        cmd = [
            "git",
            "log",
            "--format=%s",
            f"{github_before}..{github_sha}",
        ]

    # Running locally
    else:
        try:
            last_tag = subprocess.check_output(
                ["git", "describe", "--tags", "--abbrev=0"],
                text=True,
            ).strip()

            print(f"Last tag: {last_tag}")

            cmd = [
                "git",
                "log",
                "--format=%s",
                f"{last_tag}..HEAD",
            ]

        except subprocess.CalledProcessError:
            cmd = [
                "git",
                "log",
                "--format=%s",
            ]

    result = subprocess.check_output(cmd, text=True)

    return result.splitlines()


def detect_level(commits):
    level = None

    for i, msg in enumerate(commits):
        msg = msg.strip()

        print(f"{i}\t{msg}")

        if msg.startswith("MAJOR:"):
            level = "MAJOR"

        elif msg.startswith("MINOR:") and level != "MAJOR":
            level = "MINOR"

        elif msg.startswith("PATCH:") and level not in ("MAJOR", "MINOR"):
            level = "PATCH"

    return level


def bump(version, level):
    major, minor, patch = map(int, version.split("."))

    if level == "MAJOR":
        return f"{major + 1}.0.0"

    if level == "MINOR":
        return f"{major}.{minor + 1}.0"

    if level == "PATCH":
        return f"{major}.{minor}.{patch + 1}"

    return None


def main():
    content = PACKAGE_FILE.read_text()

    match = re.search(
        r'__version__\s*=\s*"([^"]+)"',
        content,
    )

    if not match:
        raise ValueError("__version__ not found")

    current = match.group(1)

    commits = get_commits()
    level = detect_level(commits)

    if level is None:
        print("No bump keyword → skip")
        return

    new_version = bump(current, level)

    updated = re.sub(
        r'(__version__\s*=\s*)"[^"]+"',
        rf'\1"{new_version}"',
        content,
    )

    PACKAGE_FILE.write_text(updated)

    print(f"{current} → {new_version} ({level})")


if __name__ == "__main__":
    main()