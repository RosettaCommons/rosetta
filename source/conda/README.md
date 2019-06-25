# Conda Build Recipes

A collection of build recipes and scripting components for building conda
packages for Rosetta components. These recipes can be used to generate
three packages:

- `pyrosetta` - A `pyrosetta` build supporting the `pyrosetta.distributed`
  namespace. Notably, this build references mandatory external
  dependencies in the PyData stack and generates specialized
  multi-threading compatible compiled components.

- `pyrosetta-binder` - A binary distribution of the `clang`-based
  `pyrosetta` binding generator used in the `pyrosetta` build process.

- `rosetta` - The rosetta suite's command-line executables and associated
  database components.

## TLDR

1. Install `docker`.
2. Use `rosetta_docker_build.sh <output dir>` to build the `rosetta` conda
   package.
3. Use `pyrosetta_docker_build.sh <output dir>` to build the `pyrosetta`
   conda package.

## Build Layout

These components are roughly organized into layers:

- `recipes`: Conda recipe definitions with (a) dependency information and
    (b) basic build scripts.

- `build`: A wrapper around conda-build rendering version information for
  a given recipe with the current tree's version information.

- `full_build`: An simple wrapper script, utilizing `version.py` to
  generate `.version.json` before invoking `build` for multiple recipes.

- `[pyrosetta|rosetta]_docker_build.sh`: Scripts using minimal build
  environments defined in linux-anvil to build broadly compatible linux
  conda packages. This should be considered the primary entrypoint to
  generate conda packages.

## Build Debugging

The `[pyrosetta|rosetta]_docker_build.sh` build scripts can be used to
diagnose failed builds within the anvil environment. To debug a failed
build:

  * Note the build prefix in the container.
    Eg.

    `/home/conda/build/rosetta_12345`

  * Re-invoke the `docker run` config used to execute the build, replacing
    the conda build call with a direct invocation of the workspace build
    script.

    Eg.
    `/home/conda/root/source/conda/build rosetta --croot /home/conda/build` 
    ->
    `bash -c 'cd /home/conda/build/rosetta_12345 && ./conda_build.sh'`
