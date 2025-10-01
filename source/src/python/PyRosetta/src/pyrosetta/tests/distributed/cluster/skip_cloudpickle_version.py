# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

import cloudpickle
import pyrosetta
import sys

from pyrosetta.utility import get_package_version


def main():
    version = get_package_version("cloudpickle")

    if not pyrosetta.rosetta.basic.was_init_called():
        pyrosetta.init(
            options="-run:constant_seed 1",
            extra_options="-out:levels core.init:0 basic.random.init_random_generator:0",
            set_logging_handler="logging",
            silent=True,
        )

    try:
        cloudpickle.dumps(pyrosetta.Pose())
        print("Confirmed that cloudpickle version {0} can pickle a Pose object.".format(version))
        sys.exit(0)
    except AttributeError as ex:
        print("Caught AttributeError in cloudpickle version {0} trying to pickle a Pose object: {1}".format(version, ex))
        sys.exit(1)

# if __name__ == "__main__":
#     print("Running: {0}".format(__file__))
#     main()
