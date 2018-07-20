import unittest
import pyrosetta
import pyrosetta.rosetta.core.pose as pose

pyrosetta.init(extra_options = "-constant_seed", set_logging_handler="logging")

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


if __name__ == "__main__":
    unittest.main()
