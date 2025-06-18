# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov
## @author Jason C. Klima

from __future__ import print_function

import glob
import json
import os
import pyrosetta
import tempfile
import unittest


class TestPyJobDistributor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
        cls.pose = pyrosetta.pose_from_sequence("DSEEKFLRRIGRFGYGYGPYE")
        cls.sf = pyrosetta.create_score_function('ref2015')
        cls.workdir = tempfile.TemporaryDirectory()
        os.chdir(cls.workdir.name)

    @classmethod
    def tearDownClass(cls):
        cls.workdir.cleanup()

    def test_jd(self):
        jd = pyrosetta.PyJobDistributor("_jd_test", 10, self.sf)

        while not jd.job_complete:
            pp = self.pose.clone()
            jd.output_decoy(pp)

        assert len( glob.glob('_jd_test*') ) == 10+1  # 10 decoys + score file

    def test_jd_at(self):
        jd = pyrosetta.PyJobDistributor("_jd_at_test", 10, self.sf, compress=True)

        while not jd.job_complete:
            pp = self.pose.clone()
            jd.output_decoy(pp)

        assert len( glob.glob('_jd_at_test*') ) == 10+1  # 10 decoys + score file

    def test_jd_serializable_scores(self):
        pdb_name = "_jd_serializable_scoretype_test"
        jd = pyrosetta.PyJobDistributor(pdb_name, 10, self.sf)

        while not jd.job_complete:
            pp = self.pose.clone()
            pp.scores["float_scoretype"] = 1e3 # Can be JSON-serialized
            pp.scores["serializable_scoretype"] = dict(foo="bar") # Can be JSON-serialized
            jd.output_decoy(pp)

        scorefile = os.path.join(self.workdir.name, "{0}.fasc".format(pdb_name))
        with open(scorefile, "r") as f:
            entries = [json.loads(line) for line in f]
        for entry in entries:
            self.assertIn("serializable_scoretype", entry, msg="PyJobDistributor did not save serializable score value.")
            self.assertEqual(entry["serializable_scoretype"]["foo"], "bar", msg="PyJobDistributor did not save serializable dictionary.")
            self.assertIn("float_scoretype", entry, msg="PyJobDistributor did not save serializable score value.")
            self.assertEqual(entry["float_scoretype"], 1e3, msg="PyJobDistributor did not save serializable float.")
            self.assertIn("total_score", entry, msg="PyJobDistributor did not save total score.")

    def test_jd_unserializable_scores(self):
        pdb_name = "_jd_unserializable_scoretype_test"
        jd = pyrosetta.PyJobDistributor(pdb_name, 10, self.sf)

        while not jd.job_complete:
            pp = self.pose.clone()
            pp.scores["unserializable_scoretype"] = 3j # Cannot be JSON-serialized
            pp.scores["serializable_scoretype"] = dict(foo="bar") # Can be JSON-serialized
            with self.assertWarns(UserWarning): # Warns about removing 'unserializable_scoretype' key
                jd.output_decoy(pp)

        scorefile = os.path.join(self.workdir.name, "{0}.fasc".format(pdb_name))
        with open(scorefile, "r") as f:
            entries = [json.loads(line) for line in f]
        for entry in entries:
            self.assertNotIn("unserializable_scoretype", entry, msg="PyJobDistributor saved unserializable score value.")
            self.assertIn("serializable_scoretype", entry, msg="PyJobDistributor did not save serializable score value.")
            self.assertEqual(entry["serializable_scoretype"]["foo"], "bar", msg="PyJobDistributor did not save serializable dictionary.")
            self.assertIn("total_score", entry, msg="PyJobDistributor did not save total score.")


if __name__ == "__main__":
    unittest.main()
