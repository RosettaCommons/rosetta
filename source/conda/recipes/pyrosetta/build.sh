#Configure!/bin/bash
#http://redsymbol.net/articles/unofficial-bash-strict-mode/

set -euo pipefail
IFS=$'\n\t'

set -x

echo "--- Env"
build_args=(
--create-package `pwd`/pyrosetta
--binder-config rosetta.config
--binder-config rosetta.distributed.config
--serialization
--multi-threaded
--no-zmq
--no-strip-module
--binder `which pyrosetta-binder`
--pybind11 `pwd`/source/external/pybind11/include
)

if test -f "`pwd`/source/.version.json"; then
  build_args+=(--version `pwd`/source/.version.json)
elif test -f "`pwd`/.release.json"; then
  build_args+=(--version `pwd`/.release.json)
fi

if [[ ! -z "${GCC:-}" ]]; then
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

  build_args+=(--compiler gcc)
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

echo "--- Build"

pushd source/src/python/PyRosetta

python build.py ${build_args[@]} -j

popd

echo "--- Setup"
pushd pyrosetta
# Run initial test to prebuild databases
python -c 'import pyrosetta; pyrosetta.init(); pyrosetta.get_score_function()(pyrosetta.pose_from_sequence("TEST"))'

pushd setup
${PREFIX}/bin/python setup.py install --single-version-externally-managed --record=record.txt > install.log
popd
