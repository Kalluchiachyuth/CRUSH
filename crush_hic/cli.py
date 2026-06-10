import os
import sys
import stat
import subprocess


def main():
    script_path = os.path.join(os.path.dirname(__file__), "CRUSH_v1.bash")

    if not os.path.isfile(script_path):
        sys.exit(f"Error: CRUSH script not found at {script_path}")

    # Ensure the bundled script is executable
    mode = os.stat(script_path).st_mode
    if not (mode & stat.S_IXUSR):
        os.chmod(script_path, mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    result = subprocess.run(["bash", script_path] + sys.argv[1:])
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
