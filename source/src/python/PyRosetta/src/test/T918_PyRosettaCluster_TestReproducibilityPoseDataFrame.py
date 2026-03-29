# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

from benchmark_test_utils import run_distributed_cluster_test_cases


run_distributed_cluster_test_cases(
    "test_reproducibility_pose_dataframe.TestReproducibilityPoseDataFrame.test_reproducibility_from_reproduce",
)
