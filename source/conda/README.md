# Conda Build Recipes

A collection of build recipes and scripting components for buliding conda
packages of Rosetta components. This can be used to generate three
packages:

- `pyrosetta` - A `pyrosetta` build supporting the `pyrosetta.distributed`
  namespace. Notably, this build references mandatory external
  dependencies in the PyData stack and generates specialized
  multi-threading compatible compiled components.

- `pyrosetta-binder` - A binary distribution of the `clang`-based
  `pyrosetta` binding generator used in the `pyrosetta` build process.

- `rosetta` - The core `rosetta_scripts`, `score`, `relax` and
  `AbinitioRelax` executables and associated database components.

## Build Layout

These components are roughly organized into layers:

- `recipes`: Conda recipe definitions with (a) dependency information and
    (b) basic build scripts.

- `build`: A wrapper around conda-build rendering version information for
  a given recipe from a pre-generated `.version.json`.

- `full_build`: An simple wrapper script, utilizing `version.py` to
  generate `.version.json` before invoking `build` for multiple recipes.

- `pipeline.yml` & `docker-compose.yml`: Definition of a standard
  "`linux-anvil`" build image, used to generate broadly compatible linux
  builds and a buildkite-based build pipeline for `linux-64` and `osx-64`
  builds deployed to a private conda channel via `rsync`.
