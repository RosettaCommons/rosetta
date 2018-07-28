# linux-anvil

This is a minimal, xenial build image based off the `conda-forge`
[`linux-anvil`](https://github.com/conda-forge/docker-images/commit/c9a9a2ae2e3f873d18212a771e7b50b53684b674).
This build image represents a best-effort attempt at a minimal build
environment for `pyrosetta` and `rosetta` builds.

The `linux-anvil` is *uniquely* designed to support `conda build` calls
while running as a non-root user. This allows image to be used to generate
build products in less-privileged contexts, writing build images into
a user-owned hosted directory via a bind mount. The target user should be
specified via the `HOST_USER_ID` env variable, rather than the standard
`-u` docker parameter.
