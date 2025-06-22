# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

import glob
import math
import os
import pyrosetta
import pyrosetta.rosetta.core.pose as pose
import tempfile
import unittest

from pyrosetta.bindings.pose import PoseScoreSerializer
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
        self.assertDictEqual(dict(test_pose.scores), {"foo" : 1})

        test_pose.scores["bar"] = 2.0
        self.assertDictEqual(dict(test_pose.scores), {"foo" : 1, "bar" : 2.0})

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
        
        # Test score value serialization
        test_pose.scores.clear()
        save_pose = pyrosetta.pose_from_sequence("SET/GET/TEST")
        save_pose.scores["test_str"] = "foo"
        save_pose.scores["test_float"] = math.pi
        save_pose.scores["test_int"] = 1111111
        type_value_dict = {
            str: "SomeString",
            float: math.pi,
            int: 42,
            dict: {"foo": 123, "bar": 456},
            tuple: (math.pi, math.pi,),
            list: ["foo", "bar", "baz"],
            set: {1, 2, 3},
            frozenset: frozenset(["a", "b", "c"]),
            bool: True,
            bytes: b"Bytes",
            bytearray: bytearray(123),
            complex: 1 + 3j,
            type(None): None,
            pyrosetta.Pose: save_pose,
            pyrosetta.rosetta.core.scoring.ScoreFunction: pyrosetta.get_score_function(),
            type(pyrosetta.rosetta.protocols.moves.NullMover): pyrosetta.rosetta.protocols.moves.NullMover, # Test `pybind11_builtins.pybind11_type` type
        }
        for obj_type, value_input in type_value_dict.items():
            # Round-trip set/get scoretype value
            if obj_type == pyrosetta.rosetta.core.scoring.ScoreFunction:
                # `ScoreFunction` instances (and some other PyRosetta instances besides `Pose()`) cannot be serialized
                with self.assertRaises(TypeError):
                    test_pose.scores[str(obj_type)] = value_input
                continue
            else:
                test_pose.scores[str(obj_type)] = value_input
                value_output = test_pose.scores[str(obj_type)]

            # Test instance types
            self.assertIsInstance(value_input, obj_type)
            if obj_type in (int, bool):
                # `int` and `bool` objects no longer cast to `float` objects during round-trip
                self.assertIsInstance(value_output, obj_type)
                self.assertNotIsInstance(value_output, float)
            else:
                self.assertIsInstance(value_output, obj_type)

            # Test values
            if obj_type == float:
                # `float` objects no longer change precision after round-trip
                self.assertEqual(value_output, value_input)
            elif obj_type == pyrosetta.Pose:
                # `Pose` objects change memory address after round-trip
                self.assertNotEqual(value_output, value_input)
                self.assertNotEqual(id(value_output), id(value_input))
                # `Pose` object properties are identical after round-trip
                self.assertEqual(value_output.size(), value_input.size()) # Test size
                self.assertEqual(value_output.annotated_sequence(), value_input.annotated_sequence()) # Test sequence
                for res in range(1, value_input.size() + 1): # Test coordinates
                    for atom in range(1, value_input.residue(res).natoms() + 1):
                        atom_input = value_input.residue(res).xyz(atom)
                        atom_output = value_output.residue(res).xyz(atom)
                        for axis in "xyz":
                            self.assertEqual(getattr(atom_input, axis), getattr(atom_output, axis))
                # `Pose` object scores are identical after round-trip
                self.assertEqual(len(value_output.scores), len(value_input.scores))
                for k in list(value_input.scores.keys()):
                    self.assertIsInstance(value_output.scores[k], type(value_input.scores[k]))
                    self.assertEqual(value_output.scores[k], value_input.scores[k])
            else:
                self.assertEqual(value_output, value_input)

        # Test custom score deserialization from SimpleMetrics data
        test_pose = test_pose.clone()
        test_pose.scores.clear()
        test_pose.scores["str"] = type_value_dict[str]
        test_pose.scores["float"] = type_value_dict[float]
        test_pose.scores["bool"] = type_value_dict[bool]
        test_pose.scores["int"] = type_value_dict[int]
        test_pose.scores["bytes"] = type_value_dict[bytes]
        test_pose.scores["complex"] = type_value_dict[complex]
        str_scores = dict(pyrosetta.rosetta.core.pose.getPoseExtraStringScores(test_pose).items())
        float_scores = dict(pyrosetta.rosetta.core.pose.getPoseExtraFloatScores(test_pose).items())
        self.assertEqual(str_scores, {}, msg="String data was not stored as SimpleMetrics data.")
        self.assertEqual(float_scores, {}, msg="Float data was not stored as SimpleMetrics data.")
        sm_data = pyrosetta.rosetta.core.simple_metrics.get_sm_data(test_pose)
        sm_real_data = dict(sm_data.get_real_metric_data())
        sm_string_data = dict(sm_data.get_string_metric_data())
        prefixes = tuple(_m.prefix for _m in PoseScoreSerializer._custom_type_metrics.values())
        self.assertIsInstance(sm_string_data["str"], str)
        self.assertFalse(sm_string_data["str"].startswith(prefixes))
        self.assertIsInstance(sm_real_data["float"], float)
        self.assertIsInstance(sm_string_data["bool"], str)
        self.assertTrue(sm_string_data["bool"].startswith(PoseScoreSerializer._custom_type_metrics["bool"].prefix))
        self.assertIsInstance(sm_string_data["int"], str)
        self.assertTrue(sm_string_data["int"].startswith(PoseScoreSerializer._custom_type_metrics["int"].prefix))
        self.assertIsInstance(sm_string_data["bytes"], str)
        self.assertTrue(sm_string_data["bytes"].startswith(PoseScoreSerializer._custom_type_metrics["bytes"].prefix))
        self.assertIsInstance(sm_string_data["complex"], str)
        self.assertTrue(sm_string_data["complex"].startswith(PoseScoreSerializer._custom_type_metrics["object"].prefix))
        # Test custom deserialization
        _custom_bool_metric = PoseScoreSerializer._custom_type_metrics["bool"]
        self.assertEqual(
            _custom_bool_metric.decode_func(sm_string_data["bool"][len(_custom_bool_metric.prefix):]),
            type_value_dict[bool],
            msg="Could not manually deserialize `bool` object."
        )
        _custom_int_metric = PoseScoreSerializer._custom_type_metrics["int"]
        self.assertEqual(
            _custom_int_metric.decode_func(sm_string_data["int"][len(_custom_int_metric.prefix):]),
            type_value_dict[int],
            msg="Could not manually deserialize `int` object."
        )
        _custom_bytes_metric = PoseScoreSerializer._custom_type_metrics["bytes"]
        self.assertEqual(
            _custom_bytes_metric.decode_func(sm_string_data["bytes"][len(_custom_bytes_metric.prefix):]),
            type_value_dict[bytes],
            msg="Could not manually deserialize `bytes` object."
        )
        _custom_arbitrary_metric = PoseScoreSerializer._custom_type_metrics["object"]
        self.assertEqual(
            _custom_arbitrary_metric.decode_func(sm_string_data["complex"][len(_custom_arbitrary_metric.prefix):]),
            type_value_dict[complex],
            msg="Could not manually deserialize `complex` object."
        )
        # Test automatic deserialization
        self.assertEqual(
            PoseScoreSerializer.maybe_decode(sm_string_data["str"]),
            type_value_dict[str],
            msg="Could not automatically deserialize `str` object."
        )
        self.assertEqual(
            PoseScoreSerializer.maybe_decode(sm_real_data["float"]),
            type_value_dict[float],
            msg="Could not automatically deserialize `float` object."
        )
        self.assertEqual(
            PoseScoreSerializer.maybe_decode(sm_string_data["bool"]),
            type_value_dict[bool],
            msg="Could not automatically deserialize `bool` object."
        )
        self.assertEqual(
            PoseScoreSerializer.maybe_decode(sm_string_data["int"]),
            type_value_dict[int],
            msg="Could not automatically deserialize `int` object."
        )
        self.assertEqual(
            PoseScoreSerializer.maybe_decode(sm_string_data["bytes"]),
            type_value_dict[bytes],
            msg="Could not automatically deserialize `bytes` object."
        )
        self.assertEqual(
            PoseScoreSerializer.maybe_decode(sm_string_data["complex"]),
            type_value_dict[complex],
            msg="Could not automatically deserialize `complex` object."
        )

        # Test deleting SimpleMetrics raising KeyError
        test_pose = test_pose.clone()
        test_pose.scores.clear()
        m = pyrosetta.rosetta.core.simple_metrics.per_residue_metrics.PerResidueClashMetric()
        m.apply(test_pose)
        with self.assertRaises(KeyError):
            test_pose.scores.pop("atomic_clashes_1")


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
