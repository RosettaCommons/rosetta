# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import argparse
import platform
import subprocess
import sys
import textwrap

from pathlib import Path


def setup_pixi_environment(env_dir):
    """
    Create a fresh pixi environment containing 'pyrosetta' and 'pyrosetta-distributed' packages.

    Note: this requires that `pixi` is an executable installed and on `${PATH}`. This function:
    - detects the current Python version
    - detects the current platform (linux/mac/windows)
    - writes a compatible 'pixi.toml' file
    - runs `pixi install` to build a new pixi environment
    """
    # Detect Python version
    py_version = f"{sys.version_info.major}.{sys.version_info.minor}"
    py_feature = f"py{sys.version_info.major}{sys.version_info.minor}"

    # Detect platform string used by GitHub Actions
    system = platform.system().lower()
    machine = platform.machine().lower()

    if system == "linux":
        plat = "linux-64" if "64" in machine else "linux-32"
    elif system == "darwin":
        plat = "osx-arm64" if "arm" in machine else "osx-64"
    elif system == "windows":
        plat = "win-64"
    else:
        raise RuntimeError(f"Unsupported platform: {system} ({machine})")

    # Build 'pixi.toml' file dynamically
    pixi_toml = textwrap.dedent(f"""
    [workspace]
    channels = ["conda-forge", "https://conda.rosettacommons.org"]
    name = "pyrosetta-pixi"
    platforms = ["{plat}"]
    version = "1.0.0"

    [dependencies]
    pyrosetta = "*"

    [pypi-dependencies]
    pyrosetta-distributed = "*"

    [feature.{py_feature}.dependencies]
    python = "{py_version}.*"

    [environments]
    {py_feature} = ["{py_feature}"]
    """)

    # Create environment directory
    env_path = Path(env_dir)
    env_path.mkdir(parents=True, exist_ok=False)
    toml_path = env_path / "pixi.toml"
    toml_path.write_text(pixi_toml.strip() + "\n")
    print(f"Created '{toml_path}' file for platform '{plat}' and Python-{py_version}")

    # Install pixi environment
    print("Running `pixi install`...")
    subprocess.run(["pixi", "install"], cwd=env_path, check=True)

    print(f"Pixi environment setup complete in directory: '{env_path}'.")


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--env_dir', type=str)
    args = parser.parse_args()
    setup_pixi_environment(args.env_dir)
