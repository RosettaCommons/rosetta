# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

from utils.distributed import exit_if_missing_numpy_requirement, run_unittest

if __name__ == "__main__":
    exit_if_missing_numpy_requirement()
    run_unittest("pyrosetta.tests.bindings.core.test_pose", timeout=60)
