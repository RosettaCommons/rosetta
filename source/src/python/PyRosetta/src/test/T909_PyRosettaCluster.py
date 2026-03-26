# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

from benchmark_test_utils import run_test_cases


run_test_cases(
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibilityMulti.test_reproducibility_from_reproduce",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibilityMulti.test_reproducibility_from_reproduce_filter_results",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibilityMulti.test_reproducibility_packer_nstruct_multi",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibilityMulti.test_reproducibility_packer_nstruct_multi_filter_results",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibilityMulti.test_reproducibility_packer_nstruct_multi_decoy_ids",
    "pyrosetta.tests.distributed.cluster.test_reproducibility.TestReproducibilityMulti.test_reproducibility_packer_nstruct_multi_decoy_ids_filter_results",
)
