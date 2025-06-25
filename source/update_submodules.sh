#!/bin/bash
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# @about This file is intended to update needed submodules prior to compilation
# To avoid issues where people don't really need all the submodules, we only update the submodules we need for
# compilation. The script optionally takes a designation of the "extras", and will conditionally update modules
# dependent on the extras settings. ( e.g. `./update_submodules.sh mpi zeromq serialization` or
# `/update_submodules.sh mpi,zeromq,serialization` )

# Work from root dir of main/ -- we have to do this for the git submodule command
cd "$( dirname "${BASH_SOURCE[0]}" )" # Magic to get to directory this file is in.
cd ..

# If we're not under git control (e.g. releases) don't bother.

if [ ! -d ".git" ]; then
  #echo "NOT UNDER GIT CONTROL"
  exit 0
fi

echo "Updating submodules needed for compilation."

# Update submodules used for all compilations
git submodule update --init -- source/external/mmtf/mmtf-cpp
git submodule update --init -- source/external/msgpack/msgpack-c.upstream
git submodule update --init -- source/external/rdkit/rdkit.upstream
git submodule update --init -- source/external/c-cmaes/c-cmaes
git submodule update --init -- source/external/gemmi_repo

# This should update all the submodules under the boost_submod directory
# On some systems Git/GitHub apparently throttles things. 
# The until loop is there to make sure this command eventually finishes.
until git submodule update --recursive --init -- source/external/boost_submod/; do
    echo ">>>>>>>>> Boost submodule update did not complete successfully. Retrying. <<<<<<<<<" | tee >(cat >&2) # Also send to stderr, for test server purposes
    sleep 10
done

# Update submodules which are extras conditional

while test $# -gt 0; do

    #echo 'Param: `'$1'`'
    # We have to consider the possibility that we're being fed a comma-separated string
    # with multiple extras. So 1) do a substring search and 2) don't stop at the first match.
    if [[ "$1" == *zeromq* ]]; then
        git submodule update --init -- source/external/libzmq
    fi
    if [[ "$1" == *tensorflow* ]]; then
        git submodule update --init -- database/protocol_data/tensorflow_graphs/tensorflow_graph_repo_submodule
    fi
    if [[ "$1" == *bcl* ]]; then
        git submodule update --init -- source/external/bcl
    fi
    #if [[ "$1" == "*FOO*" ]]; then
    #    git submodule update --init -- source/external/foo/
    #fi
    shift
done

exit 0
