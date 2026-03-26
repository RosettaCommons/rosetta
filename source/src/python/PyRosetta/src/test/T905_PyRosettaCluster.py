# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

from benchmark_test_utils import run_test_cases


run_test_cases(
    "pyrosetta.tests.distributed.cluster.test_smoke.IOTest.test_io",
    "pyrosetta.tests.distributed.cluster.test_smoke.TestInitFileSigner.test_init_file_signer",
    "pyrosetta.tests.distributed.cluster.test_smoke.PrioritiesTest.test_priorities",
    "pyrosetta.tests.distributed.cluster.test_smoke.RetriesTest.test_retries_persistent_errors",
    "pyrosetta.tests.distributed.cluster.test_smoke.RetriesTest.test_retries_succeed_on_last_retry",
    "pyrosetta.tests.distributed.cluster.test_smoke.RetriesTest.test_no_retries",
    "pyrosetta.tests.distributed.cluster.test_smoke.RetriesTest.test_retries_api",
    "pyrosetta.tests.distributed.cluster.test_smoke.TaskKeysTest.test_task_keys",
)
