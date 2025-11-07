# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import argparse
import os
import platform
import subprocess
import sys
import textwrap

from pathlib import Path


ROSETTACOMMONS_CONDA_CHANNEL = "https://conda.rosettacommons.org"


def detect_platform():
    """Detect system platform string used by GitHub Actions."""
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

    return plat


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

    # Detect platform
    plat = detect_platform()

    # Build 'pixi.toml' file dynamically
    pixi_toml = textwrap.dedent(f"""
    [workspace]
    channels = ["{ROSETTACOMMONS_CONDA_CHANNEL}", "conda-forge"]
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


def setup_uv_environment(env_dir):
    """
    Create a fresh uv environment using the 'pyrosetta-installer' package.

    Note: this requires that `uv` is an executable installed and on `${PATH}`. This function:
    - detects the current Python version
    - adds 'pyrosetta-installer' via `uv add ...`
    - runs the PyRosetta installer using `uv run python -c ...`
    """
    env_path = Path(env_dir)
    if env_path.exists():
        raise FileExistsError(f"The specified uv environment path already exists: '{env_path}'.")

    # Detect Python version
    py_version = f"{sys.version_info.major}.{sys.version_info.minor}"

    # Create uv environment using the current Python
    print(f"Creating uv environment at '{env_path}'...")
    subprocess.run(
        ["uv", "venv", str(env_path), "--python", py_version],
        check=True,
    )

    # Install pyrosetta-installer
    print("Adding 'pyrosetta-installer' and 'pip' to uv environment...")
    subprocess.run(
        ["uv", "add", "-p", str(env_path), "pyrosetta-installer"],
        check=True,
    )
    subprocess.run(
        ["uv", "add", "-p", str(env_path), "pip"],
        check=True,
    )

    # Run PyRosetta installer with mirror fallback
    print("Running PyRosetta installer in uv environment...")
    install_script = textwrap.dedent("""
        import pyrosetta_installer
        try:
            pyrosetta_installer.install_pyrosetta(
                distributed=True,
                serialization=True,
                skip_if_installed=True,
                mirror=0
            )
        except Exception as e:
            print(f"PyRosetta installation with 'mirror=0' failed: {e}. Retrying with 'mirror=1'.")
            pyrosetta_installer.install_pyrosetta(
                distributed=True,
                serialization=True,
                skip_if_installed=True,
                mirror=1
            )
    """)
    subprocess.run(
        ["uv", "run", "-p", str(env_path), "python", "-c", install_script],
        check=True,
    )

    print(f"Uv environment setup complete in directory: '{env_path}'.")


def setup_conda_environment(env_dir, env_manager="conda"):
    """
    Create a fresh conda/mamba environment containing 'pyrosetta' and 'pyrosetta-distributed' packages.

    Note: this requires that `conda` or `mamba` is an executable installed and on `${PATH}`. This function:
    - detects the current Python version
    - detects the current platform (linux/mac/windows)
    - writes a temporary 'environment.yml' file
    - runs `{env_manager} env create -f environment.yml ...` to build a new environment
    """

    # Validate environment manager
    if env_manager not in ("conda", "mamba"):
        raise ValueError("The 'env_manager' keyword argument parameter must be either 'conda' or 'mamba'.")

    # Detect Python version
    py_version = f"{sys.version_info.major}.{sys.version_info.minor}"

    # Detect platform
    plat = detect_platform()

    # Build the conda environment file dynamically
    name = os.path.basename(env_dir)
    env_yaml = textwrap.dedent(f"""
    name: {name}
    channels:
      - {ROSETTACOMMONS_CONDA_CHANNEL}
      - conda-forge
    dependencies:
      - python={py_version}.*
      - pyrosetta
      - pip
      - pip:
        - pyrosetta-distributed
    """)

    # Create environment directory
    env_path = Path(env_dir)
    env_path.mkdir(parents=True, exist_ok=False)
    env_file = env_path / "environment.yml"
    env_file.write_text(env_yaml.strip() + "\n")
    print(f"Created '{env_file}' for platform '{plat}' and Python-{py_version}")

    # Run conda/mamba to build the environment
    print(f"Running `{env_manager} env create -f {env_file} -p {env_path}`...")
    subprocess.run(
        [env_manager, "env", "create", "-f", str(env_file), "-p", str(env_path)],
        cwd=env_path,
        check=True,
    )

    print(f"{env_manager.title()} environment setup complete in directory: '{env_path}'.")


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--env_manager', type=str)
    parser.add_argument('--env_dir', type=str)
    args = parser.parse_args()
    if args.env_manager == "pixi":
        setup_pixi_environment(args.env_dir)
    elif args.env_manager == "uv":
        setup_uv_environment(args.env_dir)
    elif args.env_manager in ("conda", "mamba"):
        setup_conda_environment(args.env_dir, env_manager=args.env_manager)
    else:
        raise NotImplementedError(args.env_manager)
