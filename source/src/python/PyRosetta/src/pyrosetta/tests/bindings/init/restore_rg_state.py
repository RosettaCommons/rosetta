# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import argparse
import json
import os
import pyrosetta
import sys


def create_rg_state(init_file, json_file):
    pyrosetta.init("-run:constant_seed 1 -run:jran 1234567 -mute all", silent=True)
    rg = pyrosetta.rosetta.numeric.random.rg()
    for _ in range(999): # Cycle through MT19937 state index
        rg.uniform()
    pyrosetta.dump_init_file(init_file, verbose=False)
    value = rg.uniform() # Get uniform value
    with open(json_file, "w") as f:
        json.dump(value, f)


def restore_rg_state(init_file, json_file, restore_rg_state=None):
    assert isinstance(
        restore_rg_state, bool
    ), f"The 'restore_rg_state' keyword argument parameters must be of type `bool`. Received: {type(restore_rg_state)}"
    pyrosetta.init_from_file(
        init_file,
        restore_rg_state=restore_rg_state,
        verbose=False,
        silent=True,
    )
    rg = pyrosetta.rosetta.numeric.random.rg()
    value = rg.uniform() # Get uniform value
    with open(json_file, "r") as f:
        expected_value = json.load(f)
    if restore_rg_state:
        if value == expected_value:
            print(f"Successfully tested that RandomGenerator state was restored: {value} == {expected_value}")
            sys.exit(0) # Success
        else:
            print(f"Failed test: RandomGenerator state should be restored: {value} != {expected_value}")
            sys.exit(1) # Failure
    else:
        if value != expected_value:
            print(f"Successfully tested that RandomGenerator state was not restored: {value} != {expected_value}")
            sys.exit(0) # Success
        else:
            print(f"Failed test: RandomGenerator state should not be restored: {value} == {expected_value}")
            sys.exit(1) # Failure  


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    assert not pyrosetta.rosetta.basic.was_init_called(), "PyRosetta is already initialized!"
    parser = argparse.ArgumentParser()
    parser.add_argument('--tmp_dir', type=str)
    parser.add_argument('--create', action='store_true')
    parser.add_argument('--no-create', dest='create', action='store_false')
    parser.set_defaults(create=False)
    parser.add_argument('--restore', action='store_true')
    parser.add_argument('--no-restore', dest='restore', action='store_false')
    parser.set_defaults(restore=True)
    args = parser.parse_args()
    init_file = os.path.join(args.tmp_dir, "pyrosetta.init")
    json_file = os.path.join(args.tmp_dir, "uniform.json")
    if args.create:
        create_rg_state(init_file, json_file)
    else:
        restore_rg_state(init_file, json_file, restore_rg_state=args.restore)
