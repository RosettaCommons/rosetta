# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

from benchmark_test_utils import run_test_cases


run_test_cases(
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibility.test_reproducibility_minimizer_nstruct",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibility.test_reproducibility_minimizer_nstruct_filter_results",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibility.test_reproducibility_packer_nstruct",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibility.test_reproducibility_packer_nstruct_filter_results",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibility.test_reproducibility_packer_separate",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibility.test_reproducibility_packer_separate_filter_results",
)
