# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

from benchmark_test_utils import run_test_cases


run_test_cases(
    "pyrosetta.tests.distributed.cluster.test_smoke.WorkerPreemptionTest.test_disk_task_registry",
    "pyrosetta.tests.distributed.cluster.test_smoke.WorkerPreemptionTest.test_memory_task_registry",
    "pyrosetta.tests.distributed.cluster.test_smoke.WorkerPreemptionTest.test_disk_max_task_replicas_all",
    "pyrosetta.tests.distributed.cluster.test_smoke.WorkerPreemptionTest.test_disk_max_task_replicas_int",
)
