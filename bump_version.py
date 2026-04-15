import subprocess
import re
from pathlib import Path


def get_commits():
    try:
        last_tag = subprocess.check_output(
            ["git", "describe", "--tags", "--abbrev=0"],
            text=True
        ).strip()
        cmd = ["git", "log", f"{last_tag}..HEAD", "--pretty=%s"]
    except subprocess.CalledProcessError:
        cmd = ["git", "log", "--pretty=%s"]

    output = subprocess.check_output(cmd, text=True)
    return output.splitlines()


def detect_level(commits):
    level = None

    for c in commits:
        if c.startswith("MAJOR:"):
            return "MAJOR"
        elif c.startswith("MINOR:"):
            level = "MINOR"
        elif c.startswith("PATCH:") and level is None:
            level = "PATCH"

    return level


def bump(version, level):
    major, minor, patch = map(int, version.split("."))

    if level == "MAJOR":
        return f"{major+1}.0.0"
    elif level == "MINOR":
        return f"{major}.{minor+1}.0"
    elif level == "PATCH":
        return f"{major}.{minor}.{patch+1}"
    else:
        return None


def main():
    path = Path("oimodeler/__init__.py")
    content = path.read_text()

    match = re.search(r'__version__\s*=\s*"([^"]+)"', content)
    current = match.group(1)

    commits = get_commits()
    level = detect_level(commits)

    if level is None:
        print("No version keyword found → no bump")
        return

    new_version = bump(current, level)

    content = re.sub(
        r'(__version__\s*=\s*)"[^"]+"',
        rf'\1"{new_version}"',
        content
    )

    path.write_text(content)

    print(f"{current} → {new_version} ({level})")


if __name__ == "__main__":
    main()