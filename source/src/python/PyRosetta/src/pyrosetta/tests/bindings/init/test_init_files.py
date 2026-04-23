# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import os
import subprocess
import sys
import tempfile
import unittest

from pyrosetta.utility import has_cereal


class InitFromFileTest(unittest.TestCase):
    @unittest.skipIf(not has_cereal(), "PyRosetta is not built with serialization support.")
    def test_pipeline(self):
        """Test a PyRosetta initialization file workflow."""
        tmp_dir = tempfile.TemporaryDirectory(dir=os.getcwd(), suffix="_my_work_dir")
        cwd = os.path.dirname(__file__)
        test_scripts = [
            os.path.join(cwd, "write_test_files.py"),
            os.path.join(cwd, "dump_init_file.py"),
            os.path.join(cwd, "init_from_file.py"),
            os.path.join(cwd, "poses_from_init_file.py"),
        ]
        for test_script in test_scripts:
            cmd = "{0} {1} --tmp_dir {2}".format(sys.executable, test_script, tmp_dir.name)
            p = subprocess.run(cmd, shell=True)
            print("Return code: {0}".format(p.returncode))
            self.assertEqual(p.returncode, 0, msg=f"Test script failed: {test_script}")
        tmp_dir.cleanup()

    def test_rg_state(self):
        """Test that the RandomGenerator state can be restored using `pyrosetta.init_from_file`."""
        tmp_dir = tempfile.TemporaryDirectory(dir=os.getcwd(), suffix="_my_rg_dir")
        cwd = os.path.dirname(__file__)
        test_script = os.path.join(cwd, "restore_rg_state.py")
        base_cmd = "{0} {1} --tmp_dir {2}".format(sys.executable, test_script, tmp_dir.name)
        cmds = [
            "{0} --create".format(base_cmd),
            "{0} --no-restore".format(base_cmd),
            "{0} --restore".format(base_cmd),
        ]
        for cmd in cmds:
            p = subprocess.run(cmd, shell=True)
            print("Return code: {0}".format(p.returncode))
            self.assertEqual(p.returncode, 0, msg=f"Test script failed: {test_script}")
        tmp_dir.cleanup()


if __name__ == "__main__":
    unittest.main()
