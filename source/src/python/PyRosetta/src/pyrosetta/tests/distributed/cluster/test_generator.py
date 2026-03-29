"""
PyRosettaCluster smoke tests using the `unittest` framework.
"""
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

import os
import unittest
import uuid

from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.distributed.cluster import (
    PyRosettaCluster,
    iterate,
)
from pyrosetta.tests.distributed.cluster.unittest_utils import TestBase


class GeneratorTest(TestBase, unittest.TestCase):
    """Smoke test for the use case of the `PyRosettaCluster().generate()` method."""
    _n_tasks = 2
    _n_output_packed_poses = 2
    _parameters = (0.0, 100.0) # `float` objects for `packed_pose.update_scores` values
    _pyrosetta_kwargs = {
        "options": "-mute all",
        "extra_options": "-ex1 -multithreading:total_threads 1",
        "set_logging_handler": "logging",
        }

    @staticmethod
    def parameter_to_str(parameter):
        return f"stored_{int(parameter)}"

    @staticmethod
    def create_tasks(parameter):
        for _ in range(GeneratorTest._n_tasks):
            yield {
                "extra_options": "-ex1 -multithreading:total_threads 1",
                "set_logging_handler": "logging",
                "n_output_packed_poses": GeneratorTest._n_output_packed_poses,
                "parameter": parameter,
            }

    @staticmethod
    def my_pyrosetta_protocol_1(packed_pose, **kwargs):
        value = 1.0
        kwargs["kwargs_key_1"] = value
        packed_pose = packed_pose.update_scores(scores_key_1=value)
        yield packed_pose
        yield kwargs

    @staticmethod
    def my_pyrosetta_protocol_2(packed_pose, **kwargs):
        assert "kwargs_key_1" in kwargs
        assert "scores_key_1" in list(packed_pose.pose.scores)
        value = 2.0
        kwargs["kwargs_key_2"] = value
        parameter = kwargs.get("parameter", "string")
        packed_pose = packed_pose.update_scores(
            {
                "scores_key_2": value,
                "parameter": parameter,
                GeneratorTest.parameter_to_str(parameter): "test",
            }
        )
        packed_poses = [packed_pose.pose.clone() for _ in range(kwargs["n_output_packed_poses"])]
        return (*packed_poses, kwargs)

    @classmethod
    def get_protocols(cls):
        return [cls.my_pyrosetta_protocol_1, cls.my_pyrosetta_protocol_2]

    def test_generate_builtin_clients(self):
        """Test for `PyRosettaCluster().generate()` using built-in client instantiations."""
        instance_kwargs_update = {
            **self.instance_kwargs,
            "tasks": self.create_tasks(parameter=GeneratorTest._parameters[0]),
            "input_packed_pose": self.input_packed_pose,
            "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_builtin_clients_1_{uuid.uuid4().hex}"),
        }
        protocols = self.get_protocols()
        clients_indices = instance_kwargs_update.pop("clients_indices", None)
        resources = instance_kwargs_update.pop("resources", None)
        instance = PyRosettaCluster(**instance_kwargs_update)

        output_packed_pose_results = []
        output_kwargs_results = []
        for output_packed_pose, output_kwargs in instance.generate(
            protocols=protocols,
            clients_indices=clients_indices,
            resources=resources,
        ):
            self.assertIsInstance(output_packed_pose, PackedPose)
            self.assertIsInstance(output_kwargs, dict)
            self.assertIn("kwargs_key_1", output_kwargs)
            self.assertIn("kwargs_key_2", output_kwargs)
            self.assertIn("scores_key_1", list(output_packed_pose.pose.scores))
            self.assertIn("scores_key_2", list(output_packed_pose.pose.scores))
            instance_kwargs_update = {
                **self.instance_kwargs,
                "tasks": self.create_tasks(parameter=GeneratorTest._parameters[1]),
                "input_packed_pose": output_packed_pose,
                "protocols": protocols,
                "clients_indices": clients_indices,
                "resources": resources,
                "simulation_name": "test_generate_builtin_clients",
                "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_builtin_clients_2_{uuid.uuid4().hex}"),
            }
            for output_packed_pose, output_kwargs in iterate(**instance_kwargs_update):
                self.assertIsInstance(output_packed_pose, PackedPose)
                self.assertIsInstance(output_kwargs, dict)
                self.assertIn("kwargs_key_1", output_kwargs)
                self.assertIn("kwargs_key_2", output_kwargs)
                self.assertIn("scores_key_1", list(output_packed_pose.pose.scores))
                self.assertIn("scores_key_2", list(output_packed_pose.pose.scores))
                output_packed_pose_results.append(output_packed_pose)
                output_kwargs_results.append(output_kwargs)
        _n_results_per_parameter = (GeneratorTest._n_tasks * GeneratorTest._n_output_packed_poses)
        _n_output_results = _n_results_per_parameter ** 2
        self.assertEqual(len(output_packed_pose_results), _n_output_results)
        self.assertEqual(len(output_kwargs_results), _n_output_results)
        for packed_pose in output_packed_pose_results:
            self.assertEqual(
                packed_pose.pose.scores.get("parameter", None),
                GeneratorTest._parameters[1],
                msg="Output packed pose does not have the correct value for the 'parameter' scores key."
            )
        for kwargs in output_kwargs_results:
            self.assertEqual(
                kwargs["PyRosettaCluster_task"].get("parameter", None),
                GeneratorTest._parameters[1],
                msg="Output kwargs do not have the correct value for the 'parameter' task key."
            )
        for _parameter in GeneratorTest._parameters:
            scores_key = self.parameter_to_str(_parameter)
            self.assertEqual(
                len(list(filter(lambda packed_pose: scores_key in packed_pose.pose.scores, output_packed_pose_results))),
                _n_output_results,
                msg=f"Packed poses do not have the correct values for the '{scores_key}' scores key."
            )

    def test_generate_user_client(self):
        """
        Test for `PyRosettaCluster().generate()` with `save_all=True`,
        `dry_run=True`, and using a user-provided client.
        """
        instance_kwargs = {
            **self.instance_kwargs,
            "tasks": self.create_tasks(parameter=GeneratorTest._parameters[0]),
            "input_packed_pose": self.input_packed_pose,
            "protocols": self.get_protocols(),
            "client": self.default_client,
            "clients": None,
            "clients_indices": None,
            "resources": None,
            "save_all": True,
            "dry_run": True,
            "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_user_client_1_{uuid.uuid4().hex}"),
        }
        results = []
        for output_packed_pose, _ in iterate(**instance_kwargs):
            instance_kwargs_update = {
                **instance_kwargs,
                "input_packed_pose": output_packed_pose,
                "tasks": self.create_tasks(parameter=GeneratorTest._parameters[1]),
                "client": self.default_client, # Test passing in same client
                "simulation_name": "test_generate_user_client",
                "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_user_client_2_{uuid.uuid4().hex}"),
            }
            for result in iterate(**instance_kwargs_update):
                results.append(result)
        _n_output_packed_poses_save_all = (GeneratorTest._n_output_packed_poses + 1) # Plus one from my_pyrosetta_protocol_1
        _n_results_per_parameter = (GeneratorTest._n_tasks * _n_output_packed_poses_save_all)
        _n_output_results = _n_results_per_parameter ** 2
        self.assertEqual(
            len(results), _n_output_results, msg="Number of results with save_all failed."
        )
        self.assertListEqual(
            os.listdir(os.path.join(instance_kwargs["output_path"], instance_kwargs["decoy_dir_name"])),
            [],
            msg="Dry run failed while yielding results.",
        )
        self.assertFalse(
            os.path.isfile(os.path.join(instance_kwargs["output_path"], "scores.json")),
            msg="Dry run failed while yielding results.",
        )

    def test_generate_multi_user_clients(self):
        """
        Test for `PyRosettaCluster().generate()` with `save_all=True`,
        `dry_run=True`, and using multiple user-provided clients.
        """
        clients_indices = [0, 1]
        resources = [{"FOO": 1}, {"BAR": 9e8}]
        instance_kwargs = {
            **self.instance_kwargs,
            "tasks": self.create_tasks(parameter=GeneratorTest._parameters[0]),
            "input_packed_pose": self.input_packed_pose,
            "protocols": self.get_protocols(),
            "client": None,
            "clients": self.clients,
            "clients_indices": clients_indices,
            "resources": resources,
            "save_all": True,
            "dry_run": True,
            "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_multi_user_clients_1_{uuid.uuid4().hex}"),
        }
        results = []
        for output_packed_pose, _ in iterate(**instance_kwargs):
            instance_kwargs_update = {
                **instance_kwargs,
                "input_packed_pose": output_packed_pose,
                "tasks": self.create_tasks(parameter=GeneratorTest._parameters[1]),
                "client": None,
                "clients": self.clients, # Test passing in same clients
                "clients_indices": clients_indices,
                "resources": resources,
                "simulation_name": "test_generate_multi_user_clients",
                "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_multi_user_clients_2_{uuid.uuid4().hex}"),
            }
            for result in iterate(**instance_kwargs_update):
                results.append(result)
        _n_output_packed_poses_save_all = (GeneratorTest._n_output_packed_poses + 1) # Plus one from my_pyrosetta_protocol_1
        _n_results_per_parameter = (GeneratorTest._n_tasks * _n_output_packed_poses_save_all)
        _n_output_results = _n_results_per_parameter ** 2
        self.assertEqual(
            len(results), _n_output_results, msg="Number of results with save_all failed."
        )
        self.assertListEqual(
            os.listdir(os.path.join(instance_kwargs["output_path"], instance_kwargs["decoy_dir_name"])),
            [],
            msg="Dry run failed while yielding results.",
        )
        self.assertFalse(
            os.path.isfile(os.path.join(instance_kwargs["output_path"], "scores.json")),
            msg="Dry run failed while yielding results.",
        )

    def test_generate_partition_clients(self):
        """
        Test for `PyRosettaCluster().generate()` with `save_all=True`,
        `dry_run=True`, and partitioning user-provided clients over iterations.
        """
        resources = [{"BAZ": 1}, {"BAZ": 2}]
        instance_kwargs = {
            **self.instance_kwargs,
            "tasks": self.create_tasks(parameter=GeneratorTest._parameters[0]),
            "input_packed_pose": self.input_packed_pose,
            "protocols": self.get_protocols(),
            "client": self.clients[0], # Test passing in first client
            "clients": None,
            "clients_indices": None,
            "resources": resources,
            "save_all": True,
            "dry_run": True,
            "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_partition_clients_1_{uuid.uuid4().hex}"),
        }
        results = []
        for output_packed_pose, _ in iterate(**instance_kwargs):
            instance_kwargs_update = {
                **instance_kwargs,
                "input_packed_pose": output_packed_pose,
                "tasks": self.create_tasks(parameter=GeneratorTest._parameters[1]),
                "client": self.clients[1], # Test passing in second client
                "clients": None,
                "clients_indices": None,
                "resources": resources,
                "simulation_name": "test_generate_partition_clients",
                "output_path": os.path.join(self.workdir.name, f"outputs_test_generate_partition_clients_2_{uuid.uuid4().hex}"),
            }
            for result in iterate(**instance_kwargs_update):
                results.append(result)
        _n_output_packed_poses_save_all = (GeneratorTest._n_output_packed_poses + 1) # Plus one from my_pyrosetta_protocol_1
        _n_results_per_parameter = (GeneratorTest._n_tasks * _n_output_packed_poses_save_all)
        _n_output_results = _n_results_per_parameter ** 2
        self.assertEqual(
            len(results), _n_output_results, msg="Number of results with save_all failed."
        )
        self.assertListEqual(
            os.listdir(os.path.join(instance_kwargs["output_path"], instance_kwargs["decoy_dir_name"])),
            [],
            msg="Dry run failed while yielding results.",
        )
        self.assertFalse(
            os.path.isfile(os.path.join(instance_kwargs["output_path"], "scores.json")),
            msg="Dry run failed while yielding results.",
        )
