"""
PyRosettaCluster logger tests using the `unittest` framework.
"""
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Bobby Langan"
__email__ = "rlangan@lyell.com"

import logging
import os
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import re
import tempfile
import unittest

from pyrosetta.distributed.cluster import PyRosettaCluster


class LoggingTest(unittest.TestCase):
    _ansi_regex = re.compile(r"(?:\x1B[@-_]|[\x80-\x9F])[0-?]*[ -/]*[@-~]")

    def test_logging(self):
        """A test for capturing logging information in the distributed protocol."""

        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 300",
            set_logging_handler="logging",
        )

        def create_tasks():
            for i in range(1, 5):
                yield {
                    "extra_options": "-ex1 -multithreading:total_threads 1",
                    "seq": "LYELL" * i,
                    "set_logging_handler": "logging",
                }

        def my_pyrosetta_protocol(packed_pose, **kwargs):
            import pyrosetta
            import pyrosetta.distributed.io as io

            with self.assertLogs() as cm:
                logging.info("Logging in my_pyrosetta_protocol")
                _logger = logging.getLogger("Testing")
                _logger.info("Logging to a specific logger in my_pyrosetta_protocol")
                pyrosetta.get_score_function()

            # Testing that the logging.info syntax goes to the logger
            self.assertIn(
                "INFO:root:Logging in my_pyrosetta_protocol",
                cm.output,
                msg="root logging module not working in protocol",
            )
            # Testing that making a custom logger will rename it correctly
            self.assertIn(
                "INFO:Testing:Logging to a specific logger in my_pyrosetta_protocol",
                cm.output,
                msg="custom logging not working in protocol",
            )
            # Testing that rosetta tracers go to the logs with "set_logging_handler='logging'"
            # NOTE: This will be false if this line of the Tracer is ever changed in core.scoring.etable
            # NOTE: On some platforms there are ANSI escape codes in the Tracer
            output = [LoggingTest._ansi_regex.sub("", v) for v in cm.output]
            self.assertIn(
                "INFO:rosetta:core.scoring.etable: smooth_etable: changing atr/rep split to bottom of energy well",
                output,
                msg="rosetta logs not going to logger",
            )

            return packed_pose

        with tempfile.TemporaryDirectory() as workdir:

            output_path = os.path.join(workdir, "outputs")
            logs_dir_name = "logs"
            project_name = "PyRosettaCluster_Tests"
            simulation_name = "test"

            cluster = PyRosettaCluster(
                tasks=create_tasks,
                input_packed_pose=io.pose_from_sequence("TESTING"),
                seeds=None,
                decoy_ids=None,
                client=None,
                scheduler=None,
                scratch_dir=workdir,
                cores=None,
                processes=None,
                memory=None,
                min_workers=1,
                max_workers=1,
                nstruct=1,
                dashboard_address=None,
                compressed=True,
                logging_level="INFO",
                scorefile_name=None,
                project_name=project_name,
                simulation_name=simulation_name,
                environment=None,
                output_path=output_path,
                simulation_records_in_scorefile=False,
                decoy_dir_name="decoys",
                logs_dir_name=logs_dir_name,
                ignore_errors=False,
                timeout=0.1,
                sha1=None,
                dry_run=False,
                save_all=False,
                system_info=None,
                pyrosetta_build=None,
            )
            cluster.distribute(my_pyrosetta_protocol)

            protocol_log = os.path.join(
                output_path, logs_dir_name, project_name + "_" + simulation_name + ".log"
            )
            prc_log = os.path.join(output_path, logs_dir_name, "PyRosettaCluster.log")

            # Ensure the files populate
            if os.path.exists(protocol_log):
                self.assertGreater(
                    os.stat(protocol_log).st_size, 0, msg="Protocol log file did not populate",
                )
            else:
                print(protocol_log + " doesn't exist!")

            if os.path.exists(prc_log):
                self.assertGreater(
                    os.stat(prc_log).st_size, 0, msg="PyRosettaCluster log file did not populate",
                )
                # Ensure the PyRosettaCluster logs complete
                with open(prc_log, "r") as f:
                    last = None
                    for line in (line for line in f if line.rstrip()):
                        last = line
                    log_fields = last.split()
                    self.assertEqual(log_fields[-1], "complete!")
            else:
                print(prc_log + " doesn't exist!")


if __name__ == "__main__":
    unittest.main(verbosity=2)
