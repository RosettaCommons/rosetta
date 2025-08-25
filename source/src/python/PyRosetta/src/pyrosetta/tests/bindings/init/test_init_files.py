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


class InitFromFileTest(unittest.TestCase):
    def test_pipeline(self):
        tmp_dir = tempfile.TemporaryDirectory(dir=os.getcwd(), suffix="_my_work_dir")
        cwd = os.path.dirname(__file__)
        test_scripts = [
            os.path.join(cwd, "write_test_files.py"),
            os.path.join(cwd, "dump_init_file.py"),
            os.path.join(cwd, "init_from_file.py"),
        ]
        for test_script in test_scripts:
            cmd = "{0} {1} --tmp_dir {2}".format(sys.executable, test_script, tmp_dir.name)
            p = subprocess.run(cmd, shell=True)
            print("Return code: {0}".format(p.returncode))
            self.assertEqual(p.returncode, 0, msg=f"Test script failed: {test_script}")
        tmp_dir.cleanup()


if __name__ == "__main__":
    unittest.main()
