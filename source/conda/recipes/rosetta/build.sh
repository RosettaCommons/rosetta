#!/bin/bash

set -x
set -e

echo "--- Env"
unset MACOSX_DEPLOYMENT_TARGET

TARGET_APPS="apps"

if [[ ! -z "${GCC:-}" ]]; then
  # Override flags to just include prefix
  export CFLAGS="-I${PREFIX}/include"
  export CXXFLAGS="-I${PREFIX}/include"

  # Override flags to just include prefix
  export CFLAGS="-I${PREFIX}/include"
  export CXXFLAGS="-I${PREFIX}/include"

  # Symlink conda-provided gcc into "gcc"; pyrosetta build.py only properly
  # detects compilers named `gcc`/`g++` or `clang`/`clang++`
  mkdir -p bin
  ln -s -f ${GCC} bin/gcc
  ln -s -f ${GXX} bin/g++

  export PATH=$(pwd)/bin:$PATH
  export CC=gcc
  export CXX=g++
fi

if [[ ! -z "${CLANG:-}" ]]; then
  # override conda-provided clang compiler with system clang
  # still links against conda libc++ shared libraries
  export CLANG=/usr/bin/clang
  export CC=${CLANG}
  export CLANGXX=/usr/bin/clang++
  export CXX=${CLANGXX}
  build_args+=(--compiler ${CLANG})

  # Override flags to just include prefix
  export CFLAGS="-I${PREFIX}/include"
  export CXXFLAGS="-I${PREFIX}/include"
fi

echo "--- Configure"

pushd source/cmake
./make_project.py all

pushd build_release
cmake -G Ninja -DHDF5=ON -DCMAKE_INSTALL_PREFIX=${PREFIX}

echo "--- Build"
ninja ${TARGET_APPS}

echo "--- Install"
ninja install | grep -v Installing: | grep -v Up-to-date:

echo "--- Done"
popd
popd
