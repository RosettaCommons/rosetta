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
import pickle
import pyrosetta
import pyrosetta.rosetta.core.pose as pose
import tempfile
import unittest

from pyrosetta.bindings.scores import ClobberWarning, PoseScoreSerializer
from pyrosetta.rosetta.core.pose import (
    getPoseExtraFloatScores,
    getPoseExtraStringScores,
    setPoseExtraScore,
    hasPoseExtraScore,
    hasPoseExtraScore_str,
    clearPoseExtraScore,
    clearPoseExtraScores,
)
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta.rosetta.core.simple_metrics import TestRealMetric, TestStringMetric
from pyrosetta.rosetta.core.simple_metrics.metrics import (
    CustomRealValueMetric,
    CustomStringValueMetric,
)
from pyrosetta.rosetta.core.simple_metrics.per_residue_metrics import (
    LoadedProbabilitiesMetric,
    PerResidueClashMetric,
)
from pyrosetta.rosetta.core.simple_metrics.composite_metrics import (
    BestMutationsFromProbabilitiesMetric,
    ProtocolSettingsMetric,
)


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

    def test_score_serialization(self):
        """Test score value serialization."""
        test_pose = pyrosetta.pose_from_sequence("TESTTESTTEST")
        save_pose = pyrosetta.pose_from_sequence("SET/GET/TEST")
        save_pose.cache["test_str"] = "foo"
        save_pose.cache["test_float"] = math.pi
        save_pose.cache["test_int"] = 1111111
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
                # `ScoreFunction` instances (and some other PyRosetta instances besides `Pose` instances) cannot be serialized
                with self.assertRaises(TypeError):
                    test_pose.cache[str(obj_type)] = value_input
                continue
            else:
                test_pose.cache[str(obj_type)] = value_input
                value_output = test_pose.cache[str(obj_type)]

            # Test instance types
            self.assertIsInstance(value_input, obj_type)
            if obj_type in (int, bool):
                # `int` and `bool` objects are no longer cast to `float` objects during round-trip
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
        test_pose.cache.clear()
        test_pose.cache["str"] = type_value_dict[str]
        test_pose.cache["float"] = type_value_dict[float]
        test_pose.cache["bool"] = type_value_dict[bool]
        test_pose.cache["int"] = type_value_dict[int]
        test_pose.cache["bytes"] = type_value_dict[bytes]
        test_pose.cache["complex"] = type_value_dict[complex]
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
        test_pose.cache.clear()
        m = pyrosetta.rosetta.core.simple_metrics.per_residue_metrics.PerResidueClashMetric()
        m.apply(test_pose)
        with self.assertRaises(KeyError):
            test_pose.cache.pop("atomic_clashes_1")


class TestPoseCacheAccessor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pyrosetta.init("-ex1 -ex2aro -ex3 -run:constant_seed 1")
        cls.scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
        cls.pose = pyrosetta.pose_from_sequence("SET/GET/CACHE")
        cls.workdir = tempfile.TemporaryDirectory()
        cls.loaded_probabilities_metric = LoadedProbabilitiesMetric()
        cls.loaded_probabilities_metric.set_output_as_pdb_nums(output_as_pdb_nums=False)
        v1 = pyrosetta.rosetta.utility.vector1_std_string()
        for line in (
            "#POSNUM RESIDUETYPE WEIGHT",
            "1 ALA 0.000125",
            "1 CYS 0.0",
            "1 ASP 0.1",
            "1 GLU 0.0",
            "1 PHE 0.002197",
            "1 MET 0.972065",
            "1 ASN 0.0",
            "1 SER 0.1",
            "1 TRP 0.0",
            "2 ALA 0.0122615",
            "2 CYS 0.0",
            "2 ASP 0.0",
            "2 GLU 0.0",
            "2 PHE 0.00219777",
            "2 MET 0.9772065",
            "2 ASN 0.0",
            "2 TRP 0.0",
            "3 PRO 1.0",
            "3 THR 0.0",
            "3 TRP 0.0",
        ):
            v1.append(line)
        cls.loaded_probabilities_metric.set_probabilities_from_lines(lines=v1)
        cls.best_mutations_from_probabilities_metric = BestMutationsFromProbabilitiesMetric()
        cls.best_mutations_from_probabilities_metric.set_metric(cls.loaded_probabilities_metric)
        cls.per_residue_clash_metric = PerResidueClashMetric()
        cls.composite_settings_metric = ProtocolSettingsMetric()
        cls.custom_string_value_metric = CustomStringValueMetric()
        cls.custom_real_value_metric = CustomRealValueMetric()

    @classmethod
    def tearDownClass(cls):
        cls.workdir.cleanup()

    def test_pose_cache(self):
        self.assertEqual(self.pose.cache, {})

        # Test setting/getting energies
        self.assertEqual(self.pose.cache.energies, {})
        self.scorefxn(self.pose)
        self.assertNotEqual(self.pose.cache.energies, {})
        self.assertEqual(self.pose.cache.energies["hbond_bb_sc"], 0.0)
        self.assertEqual(self.pose.cache["hbond_sc"], 0.0)
        with self.assertRaises(NotImplementedError):
            self.pose.cache.energies.pop("total_score")
        self.pose.cache.energies.clear()
        self.assertEqual(self.pose.cache.energies, {})
        with self.assertRaises(NotImplementedError):
            self.pose.cache.energies["foo"] = 0.0

        # Test setting/getting arbitrary extra scores
        self.pose.cache.extra["foo"] = "bar"
        self.assertEqual(self.pose.cache, {"foo": "bar"})
        self.assertEqual(self.pose.cache.extra, {"foo": "bar"})
        self.assertEqual(self.pose.cache.extra["foo"], "bar")
        self.assertEqual(self.pose.cache.extra.string["foo"], "bar")
        with self.assertRaises(KeyError):
            self.pose.cache.extra.real["foo"]
        with self.assertRaises(KeyError):
            self.pose.cache.metrics["foo"]
        with self.assertRaises(KeyError):
            self.pose.cache.metrics.real["foo"]
        with self.assertRaises(KeyError):
            self.pose.cache.metrics.string["foo"]

        with self.assertWarns(UserWarning): # Float precision changes
            self.pose.cache.extra["baz"] = math.pi
        self.assertListEqual(list(self.pose.cache), ["foo", "baz"])
        self.assertListEqual(list(self.pose.cache.extra), ["foo", "baz"])
        self.assertNotEqual(self.pose.cache["baz"], math.pi, msg="Float precision was unchanged.")
        self.assertNotEqual(self.pose.cache.extra["baz"], math.pi, msg="Float precision was unchanged.")
        self.assertAlmostEqual(self.pose.cache["baz"], math.pi, places=6)
        self.assertAlmostEqual(self.pose.cache.extra["baz"], math.pi, places=6)
        self.assertAlmostEqual(self.pose.cache.extra.real["baz"], math.pi, places=6)
        self.assertEqual(self.pose.cache.extra.string["foo"], "bar")
        with self.assertRaises(KeyError):
            self.pose.cache.extra.real["foo"]
        with self.assertRaises(KeyError):
            self.pose.cache.extra.string["baz"]

        with self.assertWarns(UserWarning): # Type casting occurs
            setPoseExtraScore(self.pose, "int_to_float", 12) # Manually set
        self.assertEqual(self.pose.cache["int_to_float"], 12.0)
        self.assertEqual(self.pose.cache.extra["int_to_float"], 12.0)
        self.assertEqual(self.pose.cache.extra.real["int_to_float"], 12.0)
        with self.assertRaises(KeyError):
            self.pose.cache.extra.string["int_to_float"]

        self.pose.cache.extra["int"] = 12
        self.assertEqual(self.pose.cache["int"], 12)
        self.assertEqual(self.pose.cache.extra["int"], 12)
        self.assertEqual(self.pose.cache.extra.string["int"], 12) # Stored in string because it's serialized
        self.assertNotIsInstance(self.pose.cache["int"], float)
        with self.assertRaises(KeyError):
            self.pose.cache.extra.real["int"]

        self.pose.cache.extra["bool"] = True
        self.assertEqual(self.pose.cache["bool"], True)
        self.assertEqual(self.pose.cache.extra["bool"], True)
        self.assertEqual(self.pose.cache.extra.string["bool"], True) # Stored in string because it's serialized
        self.assertNotIsInstance(self.pose.cache["bool"], float)
        with self.assertRaises(KeyError):
            self.pose.cache.extra.real["bool"]

        setPoseExtraScore(self.pose, "bytes_to_str", b'ASCII binary') # Manually set; does not warn since it will be decoded by UTF-8
        self.assertIsInstance(self.pose.cache["bytes_to_str"], str)
        self.assertIsInstance(self.pose.cache.extra["bytes_to_str"], str)
        self.assertIsInstance(self.pose.cache.extra.string["bytes_to_str"], str)
        self.assertNotIsInstance(self.pose.cache["bytes_to_str"], bytes)
        self.assertNotIsInstance(self.pose.cache.extra["bytes_to_str"], bytes)
        self.assertNotIsInstance(self.pose.cache.extra.string["bytes_to_str"], bytes)
        self.assertEqual(self.pose.cache["bytes_to_str"], "ASCII binary")

        with self.assertWarns(UserWarning):
            setPoseExtraScore(self.pose, "invalid_bytes", pickle.dumps("raw binary")) # Manually set; warns since it's raw binary
        with self.assertRaises(UnicodeDecodeError):
            self.pose.cache["invalid_bytes"]
        with self.assertRaises(UnicodeDecodeError):
            self.pose.cache.extra["invalid_bytes"]
        with self.assertRaises(UnicodeDecodeError):
            self.pose.cache.extra.string["invalid_bytes"]
        clearPoseExtraScore(self.pose, "invalid_bytes")
        with self.assertRaises(KeyError):
            self.pose.cache["invalid_bytes"]

        for bytestring in (b'ASCII binary', pickle.dumps("raw binary")):
            self.pose.cache.extra["bytes"] = bytestring # Automatically gets serialized
            self.assertIn("bytes", self.pose.cache.extra.string)
            with self.assertRaises(KeyError):
                self.pose.cache.extra.real["bytes"]
            self.assertIsInstance(self.pose.cache["bytes"], bytes)
            self.assertEqual(self.pose.cache["bytes"], bytestring)
            self.assertEqual(self.pose.cache.extra["bytes"], bytestring)
            self.assertEqual(self.pose.cache.extra.string["bytes"], bytestring)

        self.pose.cache.extra["arbitrary"] = {1: 2, 3: b"4", 5: "6", "7": 8.0}
        self.assertIn("bytes", self.pose.cache.extra.string)
        with self.assertRaises(KeyError):
            self.pose.cache.extra.real["arbitrary"]
        self.assertDictEqual(self.pose.cache["arbitrary"], {1: 2, 3: b"4", 5: "6", "7": 8.0})
        self.assertDictEqual(self.pose.cache.extra["arbitrary"], {1: 2, 3: b"4", 5: "6", "7": 8.0})
        self.assertDictEqual(self.pose.cache.extra.string["arbitrary"], {1: 2, 3: b"4", 5: "6", "7": 8.0})

        # Test setting/getting SimpleMetric data
        self.pose.cache.metrics["my_string_metric"] = "SomeString"
        self.assertNotEqual(self.pose.cache, {"my_string_metric": "SomeString"})
        self.assertEqual(self.pose.cache.metrics, {"my_string_metric": "SomeString"})
        self.assertEqual(self.pose.cache.metrics["my_string_metric"], "SomeString")
        self.assertEqual(self.pose.cache.metrics.string["my_string_metric"], "SomeString")
        self.assertIn("my_string_metric", self.pose.cache)
        with self.assertRaises(KeyError):
            self.pose.cache.extra["my_string_metric"]
        with self.assertRaises(KeyError):
            self.pose.cache.extra.real["my_string_metric"]
        with self.assertRaises(KeyError):
            self.pose.cache.extra.string["my_string_metric"]
        with self.assertRaises(KeyError):
            self.pose.cache.metrics.real["my_string_metric"]

        self.pose.cache.metrics["my_real_metric"] = math.pi
        self.assertNotEqual(self.pose.cache, {"my_real_metric": math.pi})
        self.assertEqual(self.pose.cache.metrics, {"my_string_metric": "SomeString", "my_real_metric": math.pi})
        self.assertEqual(self.pose.cache.metrics["my_real_metric"], math.pi)
        self.assertEqual(self.pose.cache.metrics.real["my_real_metric"], math.pi)
        self.assertIn("my_real_metric", self.pose.cache)
        with self.assertRaises(KeyError):
            self.pose.cache.extra["my_real_metric"]
        with self.assertRaises(KeyError):
            self.pose.cache.extra.real["my_real_metric"]
        with self.assertRaises(KeyError):
            self.pose.cache.extra.string["my_real_metric"]
        with self.assertRaises(KeyError):
            self.pose.cache.metrics.string["my_real_metric"]

        self.custom_string_value_metric.set_value("my_data")
        self.custom_string_value_metric.apply(out_label="string_label", pose=self.pose, override_existing_data=True)
        self.assertIn("string_label", self.pose.cache)
        self.assertIn("string_label", self.pose.cache.metrics)
        self.assertIn("string_label", self.pose.cache.metrics.string)
        self.assertNotIn("string_label", self.pose.cache.metrics.real)
        self.assertNotIn("string_label", self.pose.cache.metrics.composite_real)
        self.assertNotIn("string_label", self.pose.cache.metrics.composite_string)
        self.assertNotIn("string_label", self.pose.cache.metrics.per_residue_real)
        self.assertNotIn("string_label", self.pose.cache.metrics.per_residue_string)
        self.assertNotIn("string_label", self.pose.cache.metrics.per_residue_probabilities)
        self.assertNotIn("string_label", self.pose.cache.extra)
        self.assertNotIn("string_label", self.pose.cache.extra.real)
        self.assertNotIn("string_label", self.pose.cache.extra.string)
        self.assertNotIn("string_label", self.pose.cache.energies)
        self.assertEqual(self.pose.cache["string_label"], "my_data")
        self.assertEqual(self.pose.cache.metrics["string_label"], "my_data")
        self.assertEqual(self.pose.cache.metrics.string["string_label"], "my_data")

        with self.assertRaises(TypeError): # Cannot set bytes as value in CustomRealValueMetric
            self.custom_real_value_metric.set_value(b'123')
            self.custom_real_value_metric.apply(out_label="real_label", pose=self.pose, override_existing_data=True)
        self.custom_real_value_metric.set_value(123.0)
        self.custom_real_value_metric.apply(out_label="real_label", pose=self.pose, override_existing_data=True)
        self.assertIn("real_label", self.pose.cache)
        self.assertIn("real_label", self.pose.cache.metrics)
        self.assertIn("real_label", self.pose.cache.metrics.real)
        self.assertNotIn("real_label", self.pose.cache.metrics.string)
        self.assertNotIn("real_label", self.pose.cache.metrics.composite_real)
        self.assertNotIn("real_label", self.pose.cache.metrics.composite_string)
        self.assertNotIn("real_label", self.pose.cache.metrics.per_residue_real)
        self.assertNotIn("real_label", self.pose.cache.metrics.per_residue_string)
        self.assertNotIn("real_label", self.pose.cache.metrics.per_residue_probabilities)
        self.assertNotIn("real_label", self.pose.cache.extra)
        self.assertNotIn("real_label", self.pose.cache.extra.real)
        self.assertNotIn("real_label", self.pose.cache.extra.string)
        self.assertNotIn("real_label", self.pose.cache.energies)
        self.assertEqual(self.pose.cache["real_label"], 123.0)
        self.assertEqual(self.pose.cache.metrics["real_label"], 123.0)
        self.assertEqual(self.pose.cache.metrics.real["real_label"], 123.0)

        self.custom_real_value_metric.set_value(10)
        self.custom_real_value_metric.apply(out_label="int_label", pose=self.pose, override_existing_data=True)
        self.assertNotIsInstance(self.pose.cache["int_label"], int)
        self.assertIsInstance(self.pose.cache["int_label"], float) # Integers are cast to floats when manually using CustomRealValueMetric
        self.assertEqual(self.pose.cache["int_label"], 10.0)
        self.assertIn("int_label", self.pose.cache.metrics.real) # Integers are cast to floats and stored in `pose.cache.metrics.real` when manually using CustomRealValueMetric
        self.assertNotIn("int_label", self.pose.cache.metrics.string)

        self.pose.cache["int_label"] = 10 # Saved in `pose.cache.metrics.string` because it gets serialized
        with self.assertWarns(ClobberWarning):
            self.pose.cache["int_label"] # "int_label" is in `pose.cache.metrics.real` and `pose.cache.metrics.string`
        self.pose.cache.metrics.real.pop("int_label")

        self.pose.cache["int_label"] = 10 # No clobber warning
        self.assertIsInstance(self.pose.cache["int_label"], int) # Integers are not cast to floats when using `pose.cache` setter
        self.assertNotIsInstance(self.pose.cache["int_label"], float)
        self.assertEqual(self.pose.cache["int_label"], 10)
        self.assertIn("int_label", self.pose.cache.metrics.string) # Integers are serialized into strings
        self.assertNotIn("int_label", self.pose.cache.metrics.real)

        self.composite_settings_metric.apply(out_label="settings", pose=self.pose, override_existing_data=True)
        self.assertIn("settings_constant_seed", self.pose.cache)
        self.assertIn("settings_constant_seed", self.pose.cache.metrics)
        self.assertIn("settings_constant_seed", self.pose.cache.metrics.composite_string)
        self.assertNotIn("settings_constant_seed", self.pose.cache.metrics.real)
        self.assertNotIn("settings_constant_seed", self.pose.cache.metrics.string)
        self.assertNotIn("settings_constant_seed", self.pose.cache.metrics.composite_real)
        self.assertNotIn("settings_constant_seed", self.pose.cache.metrics.per_residue_real)
        self.assertNotIn("settings_constant_seed", self.pose.cache.metrics.per_residue_string)
        self.assertNotIn("settings_constant_seed", self.pose.cache.metrics.per_residue_probabilities)
        self.assertNotIn("settings_constant_seed", self.pose.cache.extra)
        self.assertNotIn("settings_constant_seed", self.pose.cache.extra.real)
        self.assertNotIn("settings_constant_seed", self.pose.cache.extra.string)
        self.assertNotIn("settings_constant_seed", self.pose.cache.energies)
        self.assertEqual(self.pose.cache.metrics.composite_string["settings_constant_seed"], "true")
        self.assertEqual(self.pose.cache.metrics["settings_ex3"], "true")
        self.assertEqual(self.pose.cache["settings_ex2aro"], "true")

        self.loaded_probabilities_metric.apply(out_label="loaded", pose=self.pose, override_existing_data=True)
        for res in range(1, self.pose.size() + 1):
            if pyrosetta.rosetta.core.pose.res_in_chain(self.pose, resnum=res, chain="A"):
                self.assertIn(f"loaded_{res}", self.pose.cache)
                self.assertIn(f"loaded_{res}", self.pose.cache.metrics)
                self.assertIn(f"loaded_{res}", self.pose.cache.metrics.per_residue_probabilities)
                self.assertNotIn(f"loaded_{res}", self.pose.cache.metrics.real)
                self.assertNotIn(f"loaded_{res}", self.pose.cache.metrics.string)
                self.assertNotIn(f"loaded_{res}", self.pose.cache.metrics.composite_real)
                self.assertNotIn(f"loaded_{res}", self.pose.cache.metrics.composite_string)
                self.assertNotIn(f"loaded_{res}", self.pose.cache.metrics.per_residue_real)
                self.assertNotIn(f"loaded_{res}", self.pose.cache.metrics.per_residue_string)
                self.assertNotIn(f"loaded_{res}", self.pose.cache.extra)
                self.assertNotIn(f"loaded_{res}", self.pose.cache.extra.real)
                self.assertNotIn(f"loaded_{res}", self.pose.cache.extra.string)
                self.assertNotIn(f"loaded_{res}", self.pose.cache.energies)
                self.assertIsInstance(self.pose.cache[f"loaded_{res}"], dict)
                self.assertIn("TRP", self.pose.cache[f"loaded_{res}"])
                self.assertEqual(self.pose.cache[f"loaded_{res}"]["TRP"], 0.0)
                self.assertEqual(self.pose.cache.metrics[f"loaded_{res}"]["TRP"], 0.0)
                self.assertEqual(self.pose.cache.metrics.per_residue_probabilities[f"loaded_{res}"]["TRP"], 0.0)
        with self.assertRaises(RuntimeError): # Cannot override SimpleMetric with override_existing_data=False
            self.loaded_probabilities_metric.apply(out_label="loaded", pose=self.pose, override_existing_data=False)

        self.best_mutations_from_probabilities_metric.apply(out_label="best", pose=self.pose, override_existing_data=True)
        self.assertIn("best_T3P", self.pose.cache)
        self.assertIn("best_T3P", self.pose.cache.metrics)
        self.assertIn("best_T3P", self.pose.cache.metrics.composite_real)
        self.assertNotIn("best_T3P", self.pose.cache.metrics.string)
        self.assertNotIn("best_T3P", self.pose.cache.metrics.real)
        self.assertNotIn("best_T3P", self.pose.cache.metrics.composite_string)
        self.assertNotIn("best_T3P", self.pose.cache.metrics.per_residue_real)
        self.assertNotIn("best_T3P", self.pose.cache.metrics.per_residue_string)
        self.assertNotIn("best_T3P", self.pose.cache.metrics.per_residue_probabilities)
        self.assertNotIn("best_T3P", self.pose.cache.extra)
        self.assertNotIn("best_T3P", self.pose.cache.extra.real)
        self.assertNotIn("best_T3P", self.pose.cache.extra.string)
        self.assertNotIn("best_T3P", self.pose.cache.energies)
        self.assertEqual(self.pose.cache.metrics.composite_real["best_T3P"], 1.0)

        self.per_residue_clash_metric.apply(out_label="clash", pose=self.pose, override_existing_data=True)
        for res in range(1, self.pose.size() + 1):
            self.assertIn(f"clash_{res}", self.pose.cache)
            self.assertIn(f"clash_{res}", self.pose.cache.metrics)
            self.assertIn(f"clash_{res}", self.pose.cache.metrics.per_residue_real)
            self.assertNotIn(f"clash_{res}", self.pose.cache.metrics.real)
            self.assertNotIn(f"clash_{res}", self.pose.cache.metrics.string)
            self.assertNotIn(f"clash_{res}", self.pose.cache.metrics.composite_real)
            self.assertNotIn(f"clash_{res}", self.pose.cache.metrics.composite_string)
            self.assertNotIn(f"clash_{res}", self.pose.cache.metrics.per_residue_string)
            self.assertNotIn(f"clash_{res}", self.pose.cache.metrics.per_residue_probabilities)
            self.assertNotIn(f"clash_{res}", self.pose.cache.extra)
            self.assertNotIn(f"clash_{res}", self.pose.cache.extra.real)
            self.assertNotIn(f"clash_{res}", self.pose.cache.extra.string)
            self.assertNotIn(f"clash_{res}", self.pose.cache.energies)

        self.assertDictEqual(dict(self.pose.cache.metrics.per_residue_string), {})

        # Test clearing scores
        self.scorefxn(self.pose)
        ref_scores = dict(self.pose.cache.all_scores)
        self.assertDictEqual(ref_scores["metrics"]["real"], dict(self.pose.cache.metrics.real))
        self.assertDictEqual(ref_scores["metrics"]["string"], dict(self.pose.cache.metrics.string))
        self.assertDictEqual(ref_scores["metrics"]["composite_string"], dict(self.pose.cache.metrics.composite_string))
        self.assertDictEqual(ref_scores["metrics"]["composite_real"], dict(self.pose.cache.metrics.composite_real))
        self.assertDictEqual(ref_scores["metrics"]["per_residue_string"], dict(self.pose.cache.metrics.per_residue_string))
        self.assertDictEqual(ref_scores["metrics"]["per_residue_real"], dict(self.pose.cache.metrics.per_residue_real))
        self.assertDictEqual(ref_scores["metrics"]["per_residue_probabilities"], dict(self.pose.cache.metrics.per_residue_probabilities))
        self.assertDictEqual(ref_scores["extra"]["real"], dict(self.pose.cache.extra.real))
        self.assertDictEqual(ref_scores["extra"]["string"], dict(self.pose.cache.extra.string))
        self.assertDictEqual(ref_scores["energies"], dict(self.pose.cache.energies))
        ref_keys = list(self.pose.cache.all_keys)
        self.pose.cache.assert_unique_keys()

        # Test SimpleMetric data that is immutable
        with self.assertRaises(KeyError):
            self.pose.cache.metrics.composite_real.clear()
        with self.assertRaises(KeyError):
            self.pose.cache.metrics.composite_string.clear()
        with self.assertRaises(KeyError):
            self.pose.cache.metrics.per_residue_real.clear()
        with self.assertRaises(KeyError):
            self.pose.cache.metrics.per_residue_string.clear()
        with self.assertRaises(KeyError):
            self.pose.cache.metrics.per_residue_probabilities.clear()
        self.assertEqual(len(ref_keys), len(list(self.pose.cache.all_keys)))
        self.pose.cache.assert_unique_keys()
        self.pose.cache.metrics.clear() # Clear all SimpleMetric data
        self.assertGreater(len(ref_keys), len(list(self.pose.cache.all_keys)))
        self.assertEqual(len(self.pose.cache.metrics.keys()), 0)

        # Test warnings if saving objects in incorrect namespaces
        with self.assertWarns(UserWarning):
            self.pose.cache.metrics.real["str"] = "String"
        print(self.pose.cache.metrics)
        self.assertIn("str", self.pose.cache.metrics.string.keys())
        with self.assertWarns(UserWarning):
            self.pose.cache.metrics.string["float"] = 1e5
        self.assertIn("float", self.pose.cache.metrics.real.keys())
        with self.assertWarns(UserWarning):
            self.pose.cache.metrics.real["pose"] = pyrosetta.pose_from_sequence("DATA")
        self.assertIn("pose", self.pose.cache.metrics.string.keys())

        with self.assertWarns(UserWarning):
            self.pose.cache.extra.real["str"] = "String"
        with self.assertWarns(UserWarning):
            self.pose.cache.extra.string["float"] = 1e1
        with self.assertWarns(UserWarning):
            self.pose.cache.extra.real["int"] = 0
        self.assertIn("int", self.pose.cache.extra.string.keys())
        with self.assertWarns(UserWarning):
            self.pose.cache.extra.real["bytes"] = b'Bytes'
        self.assertIn("bytes", self.pose.cache.extra.string.keys())
        with self.assertWarns(UserWarning):
            self.pose.cache.extra.real["bool"] = False
        self.assertIn("bool", self.pose.cache.extra.string.keys())
        with self.assertWarns(UserWarning):
            self.pose.cache.extra.real["dict"] = dict(a=1, b=2)
        self.assertIn("dict", self.pose.cache.extra.string.keys())

        self.pose.cache.clear()
        self.pose.cache.extra.real["1"] = float(1)
        with self.assertRaises(KeyError):
            self.pose.cache.extra.string["1"] = str(1)
        self.pose.cache.extra.real.pop("1")
        self.pose.cache.extra.string["1"] = str(1)
        with self.assertRaises(KeyError):
            self.pose.cache.extra.real["1"] = float(1)

        # Test clobber warnings
        for i in range(5):
            self.pose.cache.metrics.real[str(i)] = float(i)
            self.pose.cache.metrics.string[str(i)] = str(i)
        with self.assertWarns(ClobberWarning):
            dict(self.pose.cache.metrics)

        self.pose.cache.metrics.string.clear()
        for k in self.pose.cache.metrics.real.keys():
            self.pose.cache.metrics.string[k] = "duplicate_key"
            with self.assertWarns(ClobberWarning):
                dict(self.pose.cache.metrics)
            self.pose.cache.metrics.string.pop(k)

        for k in self.pose.cache.metrics.real.keys():
            self.pose.cache.metrics.string[k] = str(self.pose.cache.metrics.real.pop(k))
        self.pose.cache.metrics.real.clear()
        for k in self.pose.cache.metrics.string.keys():
            self.pose.cache.metrics.real[k] = 1e1
            with self.assertWarns(ClobberWarning):
                dict(self.pose.cache.metrics)
            self.pose.cache.metrics.real.pop(k)

        for k in self.pose.cache.metrics.string.keys():
            self.pose.cache.metrics.real[k] = float(self.pose.cache.metrics.string[k])

        with self.assertWarns(ClobberWarning):
            dict(self.pose.cache.metrics)
        with self.assertWarns(ClobberWarning):
            dict(self.pose.cache)


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
