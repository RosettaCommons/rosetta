# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

from utils.distributed import run_distributed_cluster_test_cases


run_distributed_cluster_test_cases(
    "test_generator.GeneratorTest.test_generate_builtin_clients",
    "test_generator.GeneratorTest.test_generate_multi_user_clients",
    "test_generator.GeneratorTest.test_generate_partition_clients",
    "test_generator.GeneratorTest.test_generate_user_client",
)
