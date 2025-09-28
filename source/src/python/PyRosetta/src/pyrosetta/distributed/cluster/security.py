# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    from dask.distributed import Client, LocalCluster, Security
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.security' requires the "
        + "third-party package 'dask' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n"
    )
    raise

import logging
import os
import shutil
import subprocess

from pathlib import Path
from typing import (
    Dict,
    Generic,
    Iterable,
    NoReturn,
    Optional,
    TypeVar,
    Union,
)

from pyrosetta.distributed.cluster.hkdf import MaskedBytes, derive_task_key
from pyrosetta.distributed.cluster.serialization import NonceCache
from pyrosetta.distributed.cluster.worker_plugins import TaskSecurityPlugin


G = TypeVar("G")


class SecurityIO(Generic[G]):
    def _setup_task_security_plugin(self, clients: Dict[int, Client]) -> bytes:
        """Setup tast security worker plugin(s)."""
        prk = MaskedBytes(derive_task_key(os.urandom(32), self.instance_id))
        self._register_task_security_plugin(clients, prk)
        self.serializer.prk = prk
        assert self.serializer.__getstate__()["prk"] is None, "Pseudo-random key is not hidden on host serializer."
        if self.with_nonce:
            self.nonce_cache.prk = self.serializer.prk
            # The `NonceCache.__getstate__` method must be overridden after running `self._register_task_security_plugin()`
            # since WorkerPlugin registration first pickles `NonceCache.prk` with the default `__getstate__` method:
            self.nonce_cache.__getstate__ = NonceCache._get_state.__get__(self.nonce_cache, self.nonce_cache.__class__)
            assert self.nonce_cache.__getstate__()["prk"] is None, "Pseudo-random key is not hidden on host nonce cache."

    def _register_task_security_plugin(self, clients: Dict[int, Client], prk: MaskedBytes) -> None:
        """Register `TaskSecurityPlugin` as a dask worker plugin on dask clients."""
        for client in clients.values():
            plugin = TaskSecurityPlugin(self.instance_id, prk, self.max_nonce)
            plugin.idempotent = True # Never re-register plugin
            if hasattr(client, "register_plugin"):
                client.register_plugin(plugin=plugin, name=self.instance_id)
            else: # Deprecated since dask version 2023.9.2
                client.register_worker_plugin(plugin=plugin, name=self.instance_id, nanny=False)

    def _clients_dict_has_security(self) -> bool:
        """
        Test if the `self.clients_dict` has security enabled on all clients, excluding
        clients with `LocalCluster` clusters.
        """
        assert len(self.clients_dict) > 0, "No clients in `self.clients_dict` to test for `security` attribute."
        for _client in self.clients_dict.values():
            if not isinstance(_client.cluster, LocalCluster) and not isinstance(_client.security, Security):
                _has_security = False
                break
        else:
            _has_security = True

        return _has_security

    def _setup_with_nonce(self) -> bool:
        """Post-init hook to setup the `PyRosettaCluster().with_nonce` instance attribute."""
        if self.clients_dict:
            with_nonce = not self._clients_dict_has_security()
            if with_nonce:
                logging.warning(
                    "A PyRosettaCluster input client with a remote cluster does not have security enabled, "
                    + "so PyRosettaCluster is enabling nonce caching for security on the host and all worker "
                    + "processes. Please set the `PyRosettaCluster(max_nonce=...)` keyword argument parameter "
                    + "to limit nonce cache memory usage. To silence this warning, enable dask security on all "
                    + "PyRosettaCluster input clients that are not instances of `dask.distributed.LocalCluster`."
                )
            return with_nonce
        return not bool(self.security)
