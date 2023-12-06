# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

import os  # noqa
import pyrosetta.distributed
import pyrosetta.distributed.io as io  # noqa
import unittest

from pyrosetta.distributed.cluster import PyRosettaCluster  # noqa


class BasicTest(unittest.TestCase):
    def test_basic(self):
        def create_tasks():
            yield {
                "extra_options": "-ex1 -multithreading:total_threads 1",
                "set_logging_handler": "logging",
                "seq": "TEST",
            }

        def my_pyrosetta_protocol(packed_pose, **kwargs):
            import logging  # noqa
            import pyrosetta.distributed.io as io  # noqa

            logging.info("Testing logging from my_pyrosetta_protocol")
            # zeroDivisionError = 1 / 0  # raise intentional error

            return io.pose_from_sequence(kwargs["seq"])

        def my_pyrosetta_protocol_2(packed_pose, **kwargs):
            import logging  # noqa
            import pyrosetta.distributed.io as io  # noqa

            logging.info("Testing logging from my_pyrosetta_protocol_2")

            return packed_pose

        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )

        PyRosettaCluster(
            tasks=create_tasks,
            nstruct=1,
            seeds=[789, 123],
            decoy_ids=["0", "0"],
            scheduler=None,
            simulation_name="simulation_name",
            project_name="project_name",
            sha1=None,
            simulation_records_in_scorefile=True,
            ignore_errors=True,
        ).distribute(protocols=[my_pyrosetta_protocol, my_pyrosetta_protocol_2])


if __name__ == "__main__":
    unittest.main(verbosity=2)
