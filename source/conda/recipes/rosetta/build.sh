#!/bin/bash

set -x
set -e

unset MACOSX_DEPLOYMENT_TARGET

TARGET_APPS="rosetta_scripts score relax AbinitioRelax"

if [[ ! -z "${GCC:-}" ]]; then
  # Build via system gcc/g++ rather than conda compilers.
  # Seeing unresolvable errors on gcc 7
  export CC=`which gcc`
  export CXX=`which g++`

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

pushd source/cmake
echo "--- Build"
./make_project.py all

pushd build_release
cmake -G Ninja -DCMAKE_INSTALL_PREFIX=${PREFIX}

ninja ${TARGET_APPS}

echo "--- Install"
ninja install | grep -v Installing: | grep -v Up-to-date:
popd
popd
