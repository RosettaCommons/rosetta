#!/bin/bash
#http://redsymbol.net/articles/unofficial-bash-strict-mode/

set -euo pipefail
IFS=$'\n\t'

set -x

build_args=(
--create-package `pwd`/pyrosetta
--version `pwd`/source/.version.json
--binder-config rosetta.config
--binder-config rosetta.distributed.config
--serialization
--multi-threaded
--no-zmq
--no-strip-module
--binder `which pyrosetta-binder`
)

if [[ ! -z "${GCC:-}" ]]; then
  # Build via gcc/g++ rather than conda cc c++ compiler aliases
  # binder invokation still targets system C++ standard library
  # see linux-anvil for system gcc/g++ installation
  build_args+=(--compiler ${GCC})
  export CC=${GCC}
  export CXX=${GXX}

  # Override flags to just include prefix
  export CFLAGS="-I${PREFIX}/include"
  export CXXFLAGS="-I${PREFIX}/include"
fi

if [[ ! -z "${CLANG:-}" ]]; then
  # override conda-provided clang compiler with system clang
  # still links against conda libc++ shared libraries
  export CLANG=/usr/bin/clang
  export CC=${CLANG}
  export CLANGXX=/usr/bin/clang++
  export CXX=${CLANGXX}
  build_args+=(--compiler ${CLANG})
fi

cat source/.version.json

pushd source/src/python/PyRosetta

python build.py ${build_args[@]} -j

popd

pushd pyrosetta
# Run initial test to prebuild databases
python -c 'import pyrosetta; pyrosetta.init(); pyrosetta.get_score_function()(pyrosetta.pose_from_sequence("TEST"))'

pushd setup
${PREFIX}/bin/python setup.py install --single-version-externally-managed --record=record.txt > install.log
popd
