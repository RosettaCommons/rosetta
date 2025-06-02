# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

import glob
import os
import pyrosetta
import pyrosetta.rosetta.core.pose as pose
import tempfile
import unittest

from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta.rosetta.core.simple_metrics import TestRealMetric, TestStringMetric

pyrosetta.init(extra_options="-constant_seed", set_logging_handler="logging")

class TestPoseResidueAccessor(unittest.TestCase):

    def test_residues(self):

        pose1 = pyrosetta.pose_from_sequence('ACDEFGHI')

        # Test __len__
        self.assertEqual(0, len(pose.Pose().residues))
        self.assertEqual(8, len(pose1.residues))
        self.assertEqual(8, len(pose1))  # Deprecated

        # Test __iter__
        self.assertEqual(0, len(list(pose.Pose().residues)))
        self.assertEqual(8, len(list(pose1.residues)))
        self.assertEqual(8, len(list(pose1)))  # Deprecated

        # Test __getitem__
        # assert(pose1.residues[0] == ValueError)
        # assert(pose1.residues[0:] == ValueError)
        self.assertEqual(
            pose1.residues[1].annotated_name(), 'A[ALA:NtermProteinFull]')

        self.assertEqual(pose1.residues[6].annotated_name(), 'G')
        self.assertEqual(pose1.residues[8].annotated_name(), 'I[ILE:CtermProteinFull]')
        self.assertEqual(pose1.residues[-1].annotated_name(), 'I[ILE:CtermProteinFull]')
        self.assertEqual(pose1.residues[-3].annotated_name(), 'G')
        self.assertEqual(pose1.residues[-8].annotated_name(), 'A[ALA:NtermProteinFull]')
        self.assertEqual(
            ''.join([res.annotated_name() for res in pose.Pose().residues[:]]), '')
        self.assertEqual(
            ''.join([res.annotated_name() for res in pose1.residues[:]]),
            'A[ALA:NtermProteinFull]CDEFGHI[ILE:CtermProteinFull]')
        self.assertEqual(
            ''.join([res.annotated_name() for res in pose1.residues[1:9]]),
            'A[ALA:NtermProteinFull]CDEFGHI[ILE:CtermProteinFull]')
        self.assertEqual(
            ''.join([res.annotated_name() for res in pose1.residues[:-3]]), 'A[ALA:NtermProteinFull]CDEF')
        self.assertEqual(
            ''.join([res.annotated_name() for res in pose1.residues[3:]]), 'DEFGHI[ILE:CtermProteinFull]')
        self.assertEqual(
            ''.join([res.annotated_name() for res in pose1.residues[-6:8]]), 'DEFGH')
        self.assertEqual(
            ''.join([res.annotated_name() for res in pose1.residues[-6:8:2]]), 'DFH')
        self.assertEqual(
            ''.join([res.annotated_name() for res in pose1.residues[-6:8:3]]), 'DG')
        self.assertEqual(
            ''.join([res.annotated_name() for res in pose1[1:9]]), 'A[ALA:NtermProteinFull]CDEFGHI[ILE:CtermProteinFull]') # Deprecated
        self.assertEqual(
            ''.join([res.annotated_name() for res in pose1[-6:8]]), 'DEFGH') # Deprecated

        # #Test __iadd__
        # gly_residue = pose1.residues[6]
        # pose1.residues += gly_residue
        # self.assertEqual(
            # ''.join([res.annotated_name() for res in pose1.residues]), 'A[ALA:NtermProteinFull]CDEFGHIG[GLY:CtermProteinFull]')

        # pose2 = Pose()
        # pose2.residues += pose1.residues[1]
        # for _ in range(3):
            # pose2.residues += gly_residue
        # pose2.residues += pose1.residues[-1]
        # pose2.residues += gly_residue
        # self.assertEqual(
            # ''.join([res.annotated_name() for res in pose2.residues]), 'A[ALA:NtermProteinFull]GGGGG[GLY:CtermProteinFull]')

        # #Test __imul__
        # pose3 = Pose()
        # pose3.residues *= pose1.residues[5]
        # self.assertEqual(''.join([res.annotated_name() for res in pose3.residues]), 'F')
        # pose3.residues *= pose1.residues[5]
        # self.assertEqual(''.join([res.annotated_name() for res in pose3.residues]), 'FF')
        # pose3.residues *= pose1.residues[5]
        # self.assertEqual(''.join([res.annotated_name() for res in pose3.residues]), 'FFF')


class TestPoseScoresAccessor(unittest.TestCase):

    def test_scores(self):

        test_pose = pyrosetta.pose_from_sequence("TESTTESTTEST")

        self.assertDictEqual(dict(test_pose.scores), {})

        # Test proper overwrite of extra scores of varying types.
        test_pose.scores["foo"] = "bar"
        self.assertDictEqual(dict(test_pose.scores), {"foo" : "bar"})
        test_pose.scores["foo"] = 1
        self.assertDictEqual(dict(test_pose.scores), {"foo" : 1.0})

        test_pose.scores["bar"] = 2
        self.assertDictEqual(dict(test_pose.scores), {"foo" : 1.0, "bar" : 2.0})

        # Test score deletion
        del test_pose.scores["foo"]
        self.assertDictEqual(dict(test_pose.scores), {"bar" : 2.0})

        # Test exception when setting reserved names
        with self.assertRaises(ValueError):
            test_pose.scores["fa_atr"] = "invalid"

        # Test score update after scoring
        self.assertNotIn("fa_atr", test_pose.scores)
        pyrosetta.get_score_function()(test_pose)
        self.assertIn("fa_atr", test_pose.scores)

        # Test clear w/ energies
        with self.assertRaises(ValueError):
            del test_pose.scores["fa_atr"]

        test_pose.energies().clear()
        self.assertNotIn("fa_atr", test_pose.scores)
        self.assertDictEqual(dict(test_pose.scores), {"bar" : 2.0})

        test_pose.scores.clear()
        self.assertDictEqual(dict(test_pose.scores), dict())
        
        #Test SimpleMetric scores
        test_real = TestRealMetric()
        test_str = TestStringMetric()

        test_real.apply(test_pose)

        self.assertDictEqual(dict(test_pose.scores), {"SomeReal": 1.0})
        test_str.apply(test_pose)

        self.assertDictEqual(dict(test_pose.scores), {"SomeReal": 1.0, "SomeString": "TESTING"})

        test_pose.scores.clear()
        self.assertDictEqual(dict(test_pose.scores), dict())
        



class TestPoseResidueLabelAccessor(unittest.TestCase):

    def test_labels(self):

        test_pose = pyrosetta.pose_from_sequence("TESTTESTTEST")

        self.assertSequenceEqual(
            list(test_pose.reslabels), [set()] * len(test_pose.residues))

        test_pose.reslabels[1].add("foo")
        test_pose.reslabels[-1].add("bar")
        test_pose.reslabels[-1].add("blah")
        self.assertSequenceEqual(
            list(test_pose.reslabels),
            [{"foo"}] + [set()] * (len(test_pose.residues) - 2) + [{"bar", "blah"}])

        self.assertSequenceEqual(
            test_pose.reslabels.mask["foo"],
            [True] + [False] * (len(test_pose.residues) - 1))
        self.assertSequenceEqual(
            test_pose.reslabels.mask["bar"],
            [False] * (len(test_pose.residues) - 1) + [True])
        self.assertSequenceEqual(
            test_pose.reslabels.mask["blah"],
            [False] * (len(test_pose.residues) - 1) + [True])

        self.assertSetEqual(set(test_pose.reslabels.mask), {"foo", "bar", "blah"})

        test_pose.reslabels[-1].clear()
        self.assertSequenceEqual(
            list(test_pose.reslabels),
            [{"foo"}] + [set()] * (len(test_pose.residues) - 1))

        test_pose.reslabels[-1].add("bar")
        test_pose.reslabels[-1].add("blah")
        test_pose.reslabels[-1].discard("blah")
        self.assertSequenceEqual(
            list(test_pose.reslabels),
            [{"foo"}] + [set()] * (len(test_pose.residues) - 2) + [{"bar"}])


class TestPoseIO(unittest.TestCase):
    def test_pose_io(self):
        with tempfile.TemporaryDirectory() as workdir:
            seqs = ["TEST" * i for i in range(3, 8)]
            test_poses_seqs = [
                p.sequence() for p in pyrosetta.io.poses_from_sequences(seqs)
            ]
            self.assertListEqual(
                test_poses_seqs, seqs, msg="Sequences diverge from inputs."
            )
            self.assertListEqual(
                test_poses_seqs,
                [p.sequence() for p in pyrosetta.poses_from_sequences(tuple(seqs))],
                msg="Sequences diverge with tuple of sequences.",
            )
            for i, test_pose in enumerate(
                pyrosetta.io.poses_from_sequences(seqs), start=1
            ):
                test_pose.dump_pdb(os.path.join(workdir, "%s.pdb" % i))
            returned_poses_seqs = [
                p.sequence()
                for p in pyrosetta.poses_from_files(
                    glob.glob(os.path.join(workdir, "*.pdb"))
                )
            ]
            self.assertListEqual(
                list(sorted(test_poses_seqs)),
                list(sorted(returned_poses_seqs)),
                msg="Sequences diverge with IO.",
            )

    def test_cif_io(self):
        with tempfile.TemporaryDirectory() as workdir:
            in_pose = pyrosetta.io.pose_from_sequence("TEST/CIF")
            cif_file = os.path.join(workdir, "tmp.cif")
            pyrosetta.io.dump_cif(in_pose, cif_file)
            out_pose = pyrosetta.io.pose_from_file(cif_file)
            self.assertEqual(
                in_pose.annotated_sequence(),
                out_pose.annotated_sequence(),
                msg="Sequence recovery failed."
            )
            for res in range(1, in_pose.size() + 1):
                for atom in range(1, in_pose.residue(res).natoms() + 1):
                    atom_input = in_pose.residue(res).xyz(atom)
                    atom_output = out_pose.residue(res).xyz(atom)
                    for axis in "xyz":
                        self.assertAlmostEqual(
                            getattr(atom_input, axis),
                            getattr(atom_output, axis),
                            places=4,
                            msg="Coordinate recovery failed."
                        )

    def test_mmtf_io(self):
        with tempfile.TemporaryDirectory() as workdir:
            in_pose = pyrosetta.io.pose_from_sequence("TEST/MMTF")
            mmtf_file = os.path.join(workdir, "tmp.mmtf")
            pyrosetta.io.dump_mmtf(in_pose, mmtf_file)
            out_pose = pyrosetta.io.pose_from_file(mmtf_file)
            self.assertEqual(
                in_pose.annotated_sequence(),
                out_pose.annotated_sequence(),
                msg="Sequence recovery failed."
            )
            for res in range(1, in_pose.size() + 1):
                for atom in range(1, in_pose.residue(res).natoms() + 1):
                    atom_input = in_pose.residue(res).xyz(atom)
                    atom_output = out_pose.residue(res).xyz(atom)
                    for axis in "xyz":
                        self.assertAlmostEqual(
                            getattr(atom_input, axis),
                            getattr(atom_output, axis),
                            places=3,
                            msg="Coordinate recovery failed."
                        )

    def test_file_io(self):
        scorefxn = pyrosetta.create_score_function("ref2015")
        with tempfile.TemporaryDirectory() as workdir:
            in_pose = pyrosetta.io.pose_from_sequence("TEST/FILE")
            for filetype in ("pdb", "mmtf", "cif"):
                file = os.path.join(workdir, f"tmp.{filetype}")
                if filetype == "pdb":
                    pyrosetta.io.dump_scored_pdb(in_pose, file, scorefxn)
                else:
                    pyrosetta.io.dump_file(in_pose, file)
                out_pose = pyrosetta.io.pose_from_file(file)
                self.assertEqual(
                    in_pose.annotated_sequence(),
                    out_pose.annotated_sequence(),
                    msg="Sequence recovery failed."
                )
                for res in range(1, in_pose.size() + 1):
                    for atom in range(1, in_pose.residue(res).natoms() + 1):
                        atom_input = in_pose.residue(res).xyz(atom)
                        atom_output = out_pose.residue(res).xyz(atom)
                        for axis in "xyz":
                            self.assertAlmostEqual(
                                getattr(atom_input, axis),
                                getattr(atom_output, axis),
                                places=3,
                                msg="Coordinate recovery failed."
                            )

    def test_multimodel_pdb_io(self):
        with tempfile.TemporaryDirectory() as workdir:
            seqs = ["TESTING/" + "MANY" * i + "/SEQS" for i in range(1, 8)]
            in_poses = list(pyrosetta.io.poses_from_sequences(seqs))
            pdb_file = os.path.join(workdir, "tmp.pdb")
            pyrosetta.io.dump_multimodel_pdb(in_poses, pdb_file)
            with open(pdb_file, "r") as f:
                lines = f.readlines()
            model_lines = [line for line in lines if line.startswith("MODEL")]
            self.assertEqual(
                len(model_lines),
                len(in_poses),
                msg="Number of models recovery failed."
            )
            out_poses = list(pyrosetta.io.poses_from_multimodel_pdb(pdb_file))
            self.assertEqual(
                len(out_poses),
                len(in_poses),
                msg="Number of poses recovery failed."
            )
            for i in range(len(out_poses)):
                in_pose = in_poses[i]
                out_pose = out_poses[i]
                model = i + 1
                self.assertEqual(
                    out_pose.annotated_sequence(),
                    in_pose.annotated_sequence(),
                    msg=f"Annotated sequence recovery failed for model {model}."
                )
                for res in range(1, in_pose.size() + 1):
                    for atom in range(1, in_pose.residue(res).natoms() + 1):
                        atom_input = in_pose.residue(res).xyz(atom)
                        atom_output = out_pose.residue(res).xyz(atom)
                        for axis in "xyz":
                            self.assertAlmostEqual(
                                getattr(atom_input, axis),
                                getattr(atom_output, axis),
                                places=3,
                                msg=f"Coordinate recovery failed for model {model}."
                            )


class TestPosesToSilent(unittest.TestCase):

    def test_poses_to_silent(self):

        with tempfile.TemporaryDirectory() as workdir:

            # Tests if single pose is written to silent file and returned unchanged.
            test_pose = pyrosetta.io.pose_from_sequence("TESTTESTTEST")
            tmp_file = os.path.join(workdir, "temp.silent")
            pyrosetta.io.poses_to_silent(test_pose, tmp_file)
            returned_poses = list(pyrosetta.poses_from_silent(tmp_file))

            # Tests that the amino acid sequences does not change after 
            # being written to silent file and read back in. 
            # Cannot just assert that the two poses are equal, 
            # because writing them to the silent file adds score lines.
            self.assertEqual(
                test_pose.sequence(), 
                returned_poses[0].sequence(), 
                msg="Single sequence recovery failed.")

            # Tests that the positions of atoms are almost identical after
            # being written to silent file and read back in.
            # Rmsd will not quite be equal because the output silent file 
            # truncates xyz coordinates to the third decimal place.
            self.assertAlmostEqual(
                0., 
                all_atom_rmsd(test_pose, returned_poses[0]), 
                places=3,
                msg="Single position recovery failed.")

            # Test if a list of poses can be written to a silent file 
            # and returned unchanged.
            test_poses = [pyrosetta.io.pose_from_sequence("TEST" * i) for i in range(1, 5)]
            tmp_file = os.path.join(workdir, "temp_list.silent")
            pyrosetta.io.poses_to_silent(test_poses, tmp_file)
            returned_poses = list(pyrosetta.io.poses_from_silent(tmp_file))
            for i in range(len(test_poses)):
                self.assertEqual(
                    test_poses[i].sequence(), 
                    returned_poses[i].sequence(), 
                    msg="List sequence recovery failed.")
                self.assertAlmostEqual(
                    0.,
                    all_atom_rmsd(test_poses[i], returned_poses[i]), 
                    places=3,
                    msg="List position recovery failed.")


if __name__ == "__main__":
    unittest.main()
