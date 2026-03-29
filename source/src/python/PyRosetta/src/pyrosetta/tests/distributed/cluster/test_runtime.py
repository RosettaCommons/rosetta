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

import logging
import numpy
import sys
import time
import unittest
import uuid

from pyrosetta.distributed.cluster import PyRosettaCluster
from pyrosetta.tests.distributed.cluster.unittest_utils import (
    RuntimeTestLoggingFilter,
    TestBase,
)


class RuntimeTest(TestBase, unittest.TestCase):
    """Auxiliary runtime tests for PyRosettaCluster."""
    _pyrosetta_kwargs = {
        "options": "-mute all",
        "extra_options": "-ex1 -multithreading:total_threads 1",
        "set_logging_handler": "logging",
        }

    @staticmethod
    def create_simple_tasks(n_tasks=10):
        for i in range(n_tasks):
            yield {
                **RuntimeTest._pyrosetta_kwargs,
                "task": i,
            }

    @staticmethod
    def timing_protocol_1(packed_pose, **kwargs):
        import logging
        _logger = logging.getLogger(kwargs['PyRosettaCluster_protocol_name'])
        _logger.info(f"Running procotol {kwargs['PyRosettaCluster_protocol_name']} with task {kwargs['task']}")
        yield packed_pose

    @staticmethod
    def timing_protocol_2(packed_pose, **kwargs):
        return RuntimeTest.timing_protocol_1(packed_pose, **kwargs)

    @classmethod
    def setup_logger(cls, stream=True):
        _logger = logging.getLogger(__name__)
        if stream:
            _stream_handler = logging.StreamHandler(sys.stdout)
            _logger.addHandler(_stream_handler)
        else:
            _stream_handler = None

        _root_logger = logging.getLogger("root")
        _filter = RuntimeTestLoggingFilter()
        _root_logger.addFilter(_filter)

        _rosetta_logger = logging.getLogger("rosetta")
        _rosetta_logger.addFilter(_filter)

        _distributed_logger = logging.getLogger("pyrosetta.distributed")
        _distributed_logger.addFilter(_filter)

        return _logger, _stream_handler

    @classmethod
    def tear_down_logger(cls, _logger, _stream_handler):
        if _stream_handler is not None:
            _logger.removeHandler(_stream_handler)

    @staticmethod
    def get_mean_dt(ts):
        return numpy.mean([ts[i + 1] - ts[i] for i in range(len(ts) - 1)])

    @staticmethod
    def get_dt(t, ts):
        return t if len(ts) == 0 else (t - ts[-1])

    def test_timing_single_instance(self):
        """Runtime test with a single PyRosettaCluster instance."""
        _logger, _stream_handler = self.setup_logger(stream=True)
        _logger.info(f"Starting single PyRosettaCluster instance runtime test.")
        prc = PyRosettaCluster(
            **{
                **self.instance_kwargs,
                "input_packed_pose": self.input_packed_pose,
                "tasks": self.create_simple_tasks(),
                "clients": self.clients,
                "dry_run": True,
                "save_all": True,
                "max_delay_time": 0.0,
                "ignore_errors": False,
            }
        )
        prc_iterable = prc.generate(
            self.timing_protocol_1,
            self.timing_protocol_2,
            clients_indices=[0, 1],
            resources=[{"FOO": 1}, {"BAZ": 2}],
        )
        ts = []
        t0 = time.time()
        for i, (_, output_kwargs) in enumerate(prc_iterable):
            t = time.time() - t0
            dt = self.get_dt(t, ts)
            ts.append(t)
            task = output_kwargs.get("task")
            client_repr = output_kwargs.get("PyRosettaCluster_client_repr")
            _logger.info(f"Finished iteration {(i,)} with task {task} on client {client_repr} in {dt:0.3f} seconds")
        mean_dt = self.get_mean_dt(ts)
        _logger.info(f"Average iteration time with a single PyRosettaCluster instance: {mean_dt:0.7f} seconds")
        self.tear_down_logger(_logger, _stream_handler)

    def test_timing_multi_instance(self):
        """Runtime test for two PyRosettaCluster instances asynchronously generating results."""
        _logger, _stream_handler = self.setup_logger(stream=True)
        _logger.info(f"Starting multiple PyRosettaCluster instance runtime test.")
        setup_kwargs = {
            "dry_run": True,
            "save_all": False,
            "max_delay_time": 0.0,
            "ignore_errors": False,
        }
        prc1 = PyRosettaCluster(
            **{
                **self.instance_kwargs,
                **setup_kwargs,
                "input_packed_pose": self.input_packed_pose,
                "tasks": self.create_simple_tasks(),
                "client": self.clients[0],
            }
        )
        prc_iterate = prc1.generate(
            **{
                "protocols": self.timing_protocol_1,
                "resources": [{"BAZ": 1}],
            }
        )
        instance_kwargs_2 = {
            **self.instance_kwargs,
            **setup_kwargs,
            "client": self.clients[1],
            "simulation_name": "test_timing_multi_instance",
        }
        for k in ("protocols", "clients_indices", "resources"):
            instance_kwargs_2.pop(k, None)        
        ts = []
        t0 = time.time()
        for i, (output_packed_pose, output_kwargs) in enumerate(prc_iterate):
            t = time.time() - t0
            dt = self.get_dt(t, ts)
            ts.append(t)
            task = output_kwargs.get("task")
            client_repr = output_kwargs.get("PyRosettaCluster_client_repr")
            _logger.info(f"Finished iteration {(i,)} with task {task} on client {client_repr} in {dt:0.3f} seconds")
            prc = PyRosettaCluster(
                **{
                    **instance_kwargs_2,
                    "tasks": [{**RuntimeTest._pyrosetta_kwargs, "task": task}],
                    "input_packed_pose": output_packed_pose,
                    "simulation_name": uuid.uuid4().hex,
                }
            )
            for j, (_, output_kwargs) in enumerate(prc.generate(self.timing_protocol_2, resources=[{"BAR": 1e8}])):
                t = time.time() - t0
                dt = self.get_dt(t, ts)
                ts.append(t)
                task = output_kwargs.get("task")
                client_repr = output_kwargs.get("PyRosettaCluster_client_repr")
                _logger.info(f"Finished iteration {(i, j)} with task {task} on client {client_repr} in {dt:0.3f} seconds")
        mean_dt = self.get_mean_dt(ts)
        _logger.info(f"Average iteration time with multiple PyRosettaCluster instances: {mean_dt:0.7f} seconds")
        self.tear_down_logger(_logger, _stream_handler)
