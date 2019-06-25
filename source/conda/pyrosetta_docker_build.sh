#!/bin/bash
#http://redsymbol.net/articles/unofficial-bash-strict-mode/

set -euo pipefail
IFS=$'\n\t'

function usage() {
  echo "Usage: $BASH_SOURCE <conda-build output root>"
  exit 1
}

if [ -z "${1:-}" ]; then usage; fi

THISDIR="$( cd "$( dirname "$0" )" >/dev/null && pwd )"
ROOT=$(realpath $THISDIR/../..)
ARTIFACTS=$(realpath ${1:-})
UPLOAD_PACKAGES=False

# Conda-build copies work tree into build root, so build root must be outside
# tree to avoid recursive copy.
case $(readlink -f $(realpath $ARTIFACTS))/ 
  in $(readlink -f $ROOT)/*)
    echo "conda-build output root: '$ARTIFACTS' can not fall within working tree: '$ROOT'" && usage;
esac

# Log docker info
docker info

# In order for the conda-build process in the container to write to the mounted
# volumes, we need to run with the same id as the host machine, which is
# normally the owner of the mounted volumes, or at least has write permission
export HOST_USER_ID=$(id -u)

# Check if docker-machine is being used (normally on OSX) and get the uid from
# the VM
if hash docker-machine 2> /dev/null && docker-machine active > /dev/null; then
    export HOST_USER_ID=$(docker-machine ssh $(docker-machine active) id -u)
fi

docker build $ROOT/source/conda/linux-anvil -f $ROOT/source/conda/linux-anvil/xenial-with-gcc/Dockerfile 

DOCKER_IMAGE=$(docker build $ROOT/source/conda/linux-anvil -f $ROOT/source/conda/linux-anvil/xenial-with-gcc/Dockerfile -q)

mkdir -p "$ARTIFACTS"

set -x

docker run \
           -v "${ROOT}":/home/conda/root:rw,z \
           -v "${ARTIFACTS}":/home/conda/build:rw,z \
           -e BINSTAR_TOKEN \
           -e HOST_USER_ID=$HOST_USER_ID \
           -e UPLOAD_PACKAGES=$UPLOAD_PACKAGES \
           $DOCKER_IMAGE \
           /home/conda/root/source/conda/build pyrosetta-binder --skip-existing --croot /home/conda/build 

docker run \
           -v "${ROOT}":/home/conda/root:rw,z \
           -v "${ARTIFACTS}":/home/conda/build:rw,z \
           -e BINSTAR_TOKEN \
           -e HOST_USER_ID=$HOST_USER_ID \
           -e UPLOAD_PACKAGES=$UPLOAD_PACKAGES \
           $DOCKER_IMAGE \
           /home/conda/root/source/conda/build pyrosetta --croot /home/conda/build 
