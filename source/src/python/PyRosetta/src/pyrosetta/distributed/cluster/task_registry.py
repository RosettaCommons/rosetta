# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    import attr
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.task_registry' requires the "
        + "third-party package 'attrs' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/attrs/\n"
    )
    raise

import logging
import os
import sys

from dataclasses import asdict, dataclass
from typing import (
    Any,
    Dict,
    Generic,
    Iterator,
    Optional,
    Tuple,
    TypeVar,
    Union,
)

from pyrosetta.distributed.cluster.serialization import Serialization


G = TypeVar("G")


@dataclass
class UserArgs(Generic[G]):
    """Dataclass for the `user_spawn_thread` function argument."""
    protocol_name: str
    compressed_protocol: bytes
    compressed_packed_pose: bytes
    compressed_kwargs: bytes
    pyrosetta_init_kwargs: Dict[str, Any]
    client_repr: str
    extra_args: Dict[str, Any]
    masked_key: bytes
    task_id: str


@dataclass
class TaskRecord(Generic[G]):
    """Dataclass for PyRosettaCluster task registry entries."""
    clients_index: int
    user_args: UserArgs
    submit_kwargs: Dict[str, Any]


UnpackedTaskRecord = Tuple[int, UserArgs, Dict[str, Any]]


@attr.s(kw_only=True, slots=True, frozen=False)
class TaskRegistryBase(Generic[G]):
    """PyRosettaCluster task registry base class."""
    instance_id = attr.ib(
        type=Optional[str],
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(str)),
    )
    compression = attr.ib(
        type=Optional[Union[str, bool]],
        default="xz",
        validator=attr.validators.optional(attr.validators.instance_of((str, bool))),
    )
    serializer = attr.ib(
        type=Serialization,
        default=attr.Factory(
            lambda self: Serialization(
                instance_id=self.instance_id,
                prk=os.urandom(32),
                compression=self.compression,
                with_nonce=False,
            ),
            takes_self=True,
        ),
        validator=attr.validators.instance_of(Serialization),
        init=False,
    )

    def seal(self, task_record: TaskRecord) -> bytes:
        """Compress a task registry entry."""
        packed = self.serializer.compress_object(task_record)
        buffer = self.serializer.encoder(packed) if self.serializer.encoder else packed

        return buffer

    def unseal(self, buffer: bytes) -> TaskRecord:
        """Decompress a task registry entry."""
        packed = self.serializer.decoder(buffer) if self.serializer.decoder else buffer
        task_record = self.serializer.decompress_object(packed)

        return task_record

    def deepcopy_user_args(self, user_args: UserArgs) -> UserArgs:
        """
        Deep copy a `UserArgs` dataclass to break in-memory references to any objects
        that keep billiard subprocesses alive.
        """
        return UserArgs(**self.serializer.deepcopy_kwargs(asdict(user_args)))

    def create_task_record(
        self,
        clients_index: int,
        user_args: UserArgs,
        submit_kwargs: Dict[str, Any],
    ) -> TaskRecord:
        """Create a task record."""
        return TaskRecord(
            clients_index=clients_index,
            user_args=self.deepcopy_user_args(user_args),
            submit_kwargs=submit_kwargs,
        )

    def unpack_task_record(self, task_record: TaskRecord) -> UnpackedTaskRecord:
        """Unpack a task record."""
        return task_record.clients_index, task_record.user_args, task_record.submit_kwargs


@attr.s(kw_only=True, slots=True, frozen=False)
class DiskTaskRegistry(TaskRegistryBase[G]):
    """Task registry for on-disk PyRosettaCluster task arguments."""
    task_registry_dir = attr.ib(
        type=str,
        validator=attr.validators.instance_of(str),
    )

    def __attrs_post_init__(self) -> None:
        if not os.path.isdir(self.task_registry_dir):
            logging.info(f"Creating on-disk task registry directory: '{self.task_registry_dir}'")
            os.mkdir(self.task_registry_dir)

    def __contains__(self, key: str) -> bool:
        task_file = self._get_task_file(key, makedirs=False)

        return os.path.isfile(task_file)

    def __len__(self) -> int:
        ext = self._get_file_ext()
        total_files = 0
        for _root, _dirs, files in os.walk(self.task_registry_dir):
            total_files += sum(1 for f in files if f.endswith(ext))

        return total_files

    def __iter__(self) -> Iterator[str]:
        ext = self._get_file_ext()
        for _root, _dirs, files in os.walk(self.task_registry_dir):
            for file in files:
                if file.endswith(ext):
                    key = file[: -len(ext)]
                    yield key

    def _get_file_ext(self) -> str:
        """Get the task file extension."""
        ext = ".pkl"
        if isinstance(self.serializer.compression, str):
            ext += f".{self.serializer.compression}"

        return ext

    def _shard_key(self, key: str) -> str:
        """Shard a Dask future key for a subdirectory name."""
        key_split = key.split("-")
        if len(key_split) > 1:
            return key_split[1][:2].ljust(2, "_")
        else:
            return key[-2:].ljust(2, "_")

    def _get_task_file(self, key: str, makedirs: bool = False) -> str:
        """Get a filename for a task record."""
        cache_dir = os.path.join(self.task_registry_dir, self._shard_key(key))
        if makedirs:
            os.makedirs(cache_dir, exist_ok=True)
        ext = self._get_file_ext()

        return os.path.join(cache_dir, f"{key}{ext}")

    def total_size(self) -> int:
        """Return the total size of the on-disk task registry (in bytes)."""
        ext = self._get_file_ext()
        total_size = 0
        for root, _dirs, files in os.walk(self.task_registry_dir):
            for file in files:
                if file.endswith(ext):
                    try:
                        total_size += os.path.getsize(os.path.join(root, file))
                    except Exception:
                        continue

        return total_size

    def set(self, key: str, **kwargs: Any) -> None:
        """Set a task record into the on-disk task registry."""
        task_file = self._get_task_file(key, makedirs=True)
        if os.path.isfile(task_file):
            logging.warning(f"Task future key already exists in the on-disk task registry: '{key}'")
        with open(task_file, "wb") as f:
            f.write(self.seal(self.create_task_record(**kwargs)))

    def get(self, key: str, default: None = None) -> Optional[UnpackedTaskRecord]:
        """Get an unpacked task record from the on-disk task registry."""
        task_file = self._get_task_file(key, makedirs=False)
        if not os.path.isfile(task_file):
            logging.error(f"Task future key was not found in the on-disk task registry: '{key}'")
            return default
        with open(task_file, "rb") as f:
            try:
                task_record = self.unseal(f.read())
            except Exception as ex:
                logging.error(f"{type(ex).__name__}: Task record in the on-disk task registry is corrupted: '{key}'. {ex}")
                return default

        return self.unpack_task_record(task_record)

    def pop(self, key: str) -> None:
        """Remove a task record from the on-disk task registry."""
        task_file = self._get_task_file(key, makedirs=False)
        if os.path.isfile(task_file):
            if os.path.basename(task_file).startswith("user_spawn_thread-"):
                try:
                    os.remove(task_file)
                except Exception as ex:
                    logging.warning(f"{type(ex).__name__}: {ex}. Aborting deletion of a task registry file: '{task_file}'")
            else:
                logging.warning(
                    "Aborting deletion of a task registry file because filename does "
                    + f"not start with 'user_spawn_thread-': '{task_file}'"
                )
        else:
            logging.warning(f"Aborting deletion of a task registry file because it could not be located: '{task_file}'")

    def clear(self) -> None:
        """Clear all task records from the on-disk task registry."""
        keys = list(self)
        total_size = len(keys)
        if total_size > 0:
            logging.warning(f"Clearing {total_size} task records from the on-disk task registry: '{self.task_registry_dir}'")
            for key in keys:
                self.pop(key)
        else:
            logging.info(f"{total_size} remaining task records in the on-disk task registry.")


@attr.s(kw_only=True, slots=True, frozen=False)
class MemoryTaskRegistry(TaskRegistryBase[G]):
    """Task registry for in-memory PyRosettaCluster task arguments."""
    registry = attr.ib(
        type=Dict[str, bytes],
        default=attr.Factory(dict, takes_self=False),
        validator=attr.validators.instance_of(dict),
        init=False,
    )

    def __contains__(self, key: str) -> bool:
        return self.registry.__contains__(key)

    def __len__(self) -> int:
        return self.registry.__len__()

    def __iter__(self) -> Iterator[str]:
        return self.registry.__iter__()

    def total_size(self) -> int:
        """Return the total size of the in-memory task registry (in bytes)."""
        total_size = sys.getsizeof(self.registry)
        total_size += sum(
            sys.getsizeof(k) + sys.getsizeof(v)
            for k, v in self.registry.items()
        )

        return total_size

    def set(self, key: str, **kwargs: Any) -> None:
        """Set a task record into the in-memory task registry."""
        if key in self.registry:
            logging.warning(f"Task future key already exists in the in-memory task registry: '{key}'")
        self.registry[key] = self.seal(self.create_task_record(**kwargs))

    def get(self, key: str, default: None = None) -> Optional[UnpackedTaskRecord]:
        """Get an unpacked task record from the in-memory task registry."""
        if key not in self.registry:
            logging.error(f"Task future key was not found in the in-memory task registry: '{key}'")
            return default
        try:
            return self.unpack_task_record(self.unseal(self.registry[key]))
        except Exception as ex:
            logging.error(f"{type(ex).__name__}: Task record in the in-memory task registry is corrupted: '{key}'. {ex}")
            return default

    def pop(self, key: str) -> None:
        """Remove a task record from the in-memory task registry."""
        self.registry.pop(key, None)

    def clear(self) -> None:
        """Clear all task records from the in-memory task registry."""
        total_size = len(self)
        if total_size > 0:
            logging.warning(f"Clearing {total_size} task records from the in-memory task registry.")
            self.registry.clear()
        else:
            logging.info(f"{total_size} remaining task records in the in-memory task registry.")
