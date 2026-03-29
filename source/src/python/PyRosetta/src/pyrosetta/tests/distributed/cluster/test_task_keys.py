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

import inspect
import os
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import sys
import tempfile
import unittest
import uuid

from pyrosetta.distributed.cluster import run
from pyrosetta.distributed.cluster.exceptions import WorkerError
from pyrosetta.distributed.cluster.multiprocessing import run_protocol, user_protocol


@unittest.skipIf(
    sys.version_info[:2] < (3, 8),
    "Unit test implements positional-only argument separators in function signatures, which is supported by Python-3.8+ only."
)
class TaskKeysTest(unittest.TestCase):
    """Smoke tests for user-defined task dictionary keys and user-defined PyRosetta protocol input signatures."""
    @classmethod
    def setUpClass(cls):
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 200",
            set_logging_handler="logging",
        )
        cls.input_packed_pose = io.pose_from_sequence("ARG/SEP")
        cls.workdir = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.workdir.cleanup()

    @staticmethod
    def positional_arbitrary_parameter_protocol(baz, /, **kwargs):
        """Test Case #1: A user-defined PyRosetta protocol with a positional-only parameter named 'baz'."""
        return baz, kwargs

    @staticmethod
    def positional_or_keyword_arbitrary_parameter_protocol(baz, **kwargs):
        """Test Case #2: A user-defined PyRosetta protocol with a positional-or-keyword parameter named 'baz'."""
        return baz, kwargs

    @staticmethod
    def positional_custom_parameter_protocol(foo, /, **kwargs):
        """Test Case #3: A user-defined PyRosetta protocol with a positional-only parameter named 'foo'."""
        return foo, kwargs

    @staticmethod
    def positional_or_keyword_custom_parameter_protocol(foo, **kwargs):
        """Test Case #4: A user-defined PyRosetta protocol with a positional-or-keyword parameter named 'foo'."""
        return foo, kwargs

    @staticmethod
    def positional_parameter_protocol(packed_pose, /, **kwargs):
        """Test Case #5: A user-defined PyRosetta protocol with a positional-only parameter named 'packed_pose'."""
        return packed_pose, kwargs

    @staticmethod
    def positional_or_keyword_parameter_protocol(packed_pose, **kwargs):
        """Test Case #6: A user-defined PyRosetta protocol with a positional-or-keyword parameter named 'packed_pose'."""
        return packed_pose, kwargs

    def test_task_keys(self):
        """Test compatibility between user-defined task dictionary keys and user-defined PyRosetta protocol input signatures."""
        task_keys = set()
        for func in (run_protocol, user_protocol):
            sig = inspect.signature(func)
            for name, param in sig.parameters.items():
                if str(param.kind) != "VAR_KEYWORD": # Add all except variadic keyword parameters
                    task_keys.add(name)
        task_keys.update({"foo", "bar"}) # Add custom task keys
        instance_kwargs = dict(
            tasks=dict(map(reversed, enumerate(task_keys))),
            input_packed_pose=self.input_packed_pose,
            scheduler=None,
            cores=None,
            processes=None,
            memory=None,
            min_workers=1,
            max_workers=1,
            nstruct=1,
            project_name="PyRosettaCluster_Tests",
            simulation_name=uuid.uuid4().hex,
            scratch_dir=os.path.join(self.workdir.name, "scratch"),
            sha1=None,
        )
        output_path = os.path.join(self.workdir.name, "test_task_keys")
        _sep = "*" * 60

        # --- Test case #1 ---
        # Test a user-defined PyRosetta protocol defining a positional-only parameter named 'baz':
        # The user-defined task dictionary does not contain the 'baz' key, and the user-defined PyRosetta
        # protocol defines a positional-only parameter named 'baz', so upon calling the user-defined PyRosetta protocol,
        # PyRosettaCluster passes the `packed_pose` variable as a positional argument to the 'baz' parameter, and the
        # 'packed_pose' keyword argument from the unpacked user-defined task dictionary is collected by the variadic
        # keyword parameter, which avoids a `TypeError` since it's valid Python syntax.
        protocol = TaskKeysTest.positional_arbitrary_parameter_protocol
        parameter_name = "baz"
        self.assertNotIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_1")

        # --- Test case #2 ---
        # Test a user-defined PyRosetta protocol defining a positional-or-keyword parameter named 'baz':
        # The user-defined task dictionary does not contain the 'baz' key, and the user-defined PyRosetta
        # protocol defines a positional-or-keyword parameter named 'baz', so upon calling the user-defined PyRosetta protocol,
        # PyRosettaCluster passes the `packed_pose` variable as a positional argument to the 'baz' parameter, and the
        # 'packed_pose' keyword argument from the unpacked user-defined task dictionary is collected by the variadic
        # keyword parameter, which avoids a `TypeError` since it's valid Python syntax.
        protocol = TaskKeysTest.positional_or_keyword_arbitrary_parameter_protocol
        parameter_name = "baz"
        self.assertNotIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_2")

        # --- Test case #3 ---
        # Test a user-defined PyRosetta protocol defining a positional-only parameter named 'foo':
        # Although the user-defined task dictionary contains the 'foo' key, because the user-defined PyRosetta
        # protocol defines a positional-only parameter named 'foo', upon calling the user-defined PyRosetta protocol,
        # PyRosettaCluster passes the `packed_pose` variable as a positional argument to the 'foo' parameter, and the
        # 'foo' keyword argument from the unpacked user-defined task dictionary is collected by the variadic
        # keyword parameter, which avoids a `TypeError` since it's valid Python syntax.
        protocol = TaskKeysTest.positional_custom_parameter_protocol
        parameter_name = "foo"
        self.assertIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_3")

        # --- Test case #4 ---
        # Test a user-defined PyRosetta protocol defining a positional-or-keyword parameter named 'foo':
        protocol = TaskKeysTest.positional_or_keyword_custom_parameter_protocol
        parameter_name = "foo"
        self.assertIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        print(f"{_sep} Begin testing expected `TypeError` in billiard subprocess {_sep}", flush=True)
        with self.assertRaises(WorkerError):
            # Catch exception in billiard subprocess raised due to the user-defined task dictionary containing the
            # 'foo' key and the user-defined PyRosetta protocol defining the positional-or-keyword parameter 'foo'.
            # When calling the user-defined PyRosetta protocol, PyRosettaCluster passes the `packed_pose` variable as
            # a positional argument to the 'foo' parameter, and the user-defined task dictionary is unpacked with the
            # 'foo' keyword argument name, causing the exception:
            #     TypeError: TaskKeysTest.positional_or_keyword_custom_parameter_protocol() got multiple values for argument 'foo'
            # However, this is a user error. To fix it, the user must either: (i) remove the 'foo' key from the user-defined
            # task dictionary; (ii) add a positional-only argument separator to the user-defined PyRosetta protocol input
            # signature to make the 'foo' parameter positional-only (as in Test Case #3); or (iii) change the positional-or-keyword
            # parameter name in the user-defined PyRosetta protocol input signature so that it doesn't clash with the 'foo' key
            # from the user-defined task dictionary (as in Test Case #2).
            run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_4")
        print(f"{_sep} End testing expected `TypeError` in billiard subprocess {_sep}", flush=True)

        # --- Test case #5 ---
        # Test a user-defined PyRosetta protocol defining a positional-only parameter named 'packed_pose':
        # Although the user-defined task dictionary contains the 'packed_pose' key, because the user-defined PyRosetta
        # protocol defines a positional-only parameter named 'packed_pose', upon calling the user-defined PyRosetta protocol,
        # PyRosettaCluster passes the `packed_pose` variable as a positional argument to the 'packed_pose' parameter, and the
        # 'packed_pose' keyword argument from the unpacked user-defined task dictionary is collected by the variadic
        # keyword parameter, which avoids a `TypeError` since it's valid Python syntax.
        protocol = TaskKeysTest.positional_parameter_protocol
        parameter_name = "packed_pose"
        self.assertIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_5")

        # --- Test case #6 ---
        # Test a user-defined PyRosetta protocol defining a positional-or-keyword parameter named 'packed_pose':
        protocol = TaskKeysTest.positional_or_keyword_parameter_protocol
        parameter_name = "packed_pose"
        self.assertIn(parameter_name, instance_kwargs.get("tasks"))
        self.assertIn(parameter_name, inspect.signature(protocol).parameters)
        print(f"{_sep} Begin testing expected `TypeError` in billiard subprocess {_sep}", flush=True)
        with self.assertRaises(WorkerError):
            # Catch exception in billiard subprocess raised due to the user-defined task dictionary containing the
            # 'packed_pose' key and the user-defined PyRosetta protocol defining the positional-or-keyword parameter 'packed_pose'.
            # When calling the user-defined PyRosetta protocol, PyRosettaCluster passes the `packed_pose` variable as a positional
            # argument to the 'packed_pose' parameter, and the user-defined task dictionary is unpacked with the 'packed_pose'
            # keyword argument name, causing the exception:
            #     TypeError: TaskKeysTest.positional_or_keyword_parameter_protocol() got multiple values for argument 'packed_pose'
            # However, this is a user error. To fix it, the user must either: (i) remove the 'packed_pose' key from the
            # user-defined task dictionary; (ii) add a positional-only argument separator to the user-defined PyRosetta protocol
            # input signature to make the 'packed_pose' parameter positional-only (as in Test Case #5); or (iii) change the
            # positional-or-keyword parameter name in the user-defined PyRosetta protocol input signature so that it doesn't clash
            # with the 'packed_pose' key from the user-defined task dictionary (as in Test Case #2).
            run(**instance_kwargs, protocols=protocol, output_path=f"{output_path}_6")
        print(f"{_sep} End testing expected `TypeError` in billiard subprocess {_sep}", flush=True)
