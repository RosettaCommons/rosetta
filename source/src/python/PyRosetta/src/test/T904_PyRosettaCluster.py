# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

from benchmark_test_utils import run_test_cases


run_test_cases(
    "pyrosetta.tests.distributed.cluster.test_smoke.MultipleClientsTest.test_clients",
    "pyrosetta.tests.distributed.cluster.test_smoke.ResourcesTest.test_resources",
    "pyrosetta.tests.distributed.cluster.test_smoke.ResourcesTest.test_resources_clients",
    "pyrosetta.tests.distributed.cluster.test_smoke.GeneratorTest.test_generate_builtin_clients",
    "pyrosetta.tests.distributed.cluster.test_smoke.GeneratorTest.test_generate_multi_user_clients",
    "pyrosetta.tests.distributed.cluster.test_smoke.GeneratorTest.test_generate_partition_clients",
    "pyrosetta.tests.distributed.cluster.test_smoke.GeneratorTest.test_generate_user_client",
    "pyrosetta.tests.distributed.cluster.test_smoke.RuntimeTest.test_timing_multi_instance",
    "pyrosetta.tests.distributed.cluster.test_smoke.RuntimeTest.test_timing_single_instance",
)
