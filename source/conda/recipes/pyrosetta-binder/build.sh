#!/bin/bash
#http://redsymbol.net/articles/unofficial-bash-strict-mode/

set -exo pipefail
IFS=$'\n\t'

if [[ ! -z "${CLANG:-}" ]]; then
  # override conda-provided clang compiler with system clang
  # still links against conda libc++ shared libraries
  export CLANG=/usr/bin/clang
  export CC=${CLANG}
  export CLANGXX=/usr/bin/clang++
  export CXX=${CLANGXX}
fi

python build.py -j

LIB="${PREFIX}/lib/pyrosetta-binder.${PKG_VERSION}"

mkdir -p $LIB
cp build/llvm-4.0.0/build_*/bin/binder ${LIB}

mkdir -p ${LIB}/clang/lib
cp -r build/llvm-4.0.0/tools/clang/lib/Headers ${LIB}/clang/lib/Headers

mkdir -p ${LIB}/pybind11/
cp -r build/pybind11/include ${LIB}/pybind11/include

mkdir -p ${PREFIX}/bin
pushd ${PREFIX}/bin
cat > pyrosetta-binder.${PKG_VERSION} <<SHIM
#!/usr/bin/env python3

import sys
import os
import subprocess

real_bin = os.path.realpath(__file__)
prefix = os.path.dirname(os.path.dirname(real_bin))
binder_install_path = os.path.join(prefix, "lib", os.path.basename(real_bin))

if len(sys.argv) > 1 and sys.argv[1] == "--pybind11-include-path":
    assert len(sys.argv) == 2
    print(os.path.join(binder_install_path, "pybind11/include"))
else:
    if not "--" in sys.argv:
        sys.argv.append("--")

    subprocess.call(
        [os.path.join(binder_install_path, "binder")] +
        sys.argv[1:] + 
        ["-I" + os.path.join(binder_install_path, "clang/lib/Headers")]
    )
SHIM

chmod a+x pyrosetta-binder.${PKG_VERSION}

ln -f -s pyrosetta-binder.${PKG_VERSION} pyrosetta-binder
