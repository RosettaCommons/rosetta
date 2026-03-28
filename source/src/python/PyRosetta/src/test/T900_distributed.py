# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

from benchmark_test_utils import run_test_cases


run_test_cases(
    "pyrosetta.tests.bindings.init.test_init_files",
    "pyrosetta.tests.bindings.core.test_pose",
    "pyrosetta.tests.distributed.test_concurrency",
    "pyrosetta.tests.distributed.test_dask",
    "pyrosetta.tests.distributed.test_gil",
    "pyrosetta.tests.distributed.test_smoke",
    "pyrosetta.tests.distributed.test_viewer",
    "pyrosetta.tests.numeric.test_alignment",
)
