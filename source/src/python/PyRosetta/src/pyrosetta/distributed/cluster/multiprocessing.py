# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import billiard
    from dask.distributed import get_client, get_worker
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.multiprocessing' requires the "
        + "third-party packages 'billiard' and 'dask.distributed' as dependencies!\n"
        + "Please install the packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/billiard/\n"
        + "https://pypi.org/project/distributed/\n"
    )
    raise

import tempfile
import time

from pyrosetta.distributed import requires_init
from pyrosetta.distributed.cluster.base import (
    _get_residue_type_set,
    capture_task_metadata,
    TASK_ARGS_PLUGIN_NAME,
)
from pyrosetta.distributed.cluster.converters import (
    _parse_empty_queue,
    _parse_protocol_results,
    _parse_target_results,
)
from pyrosetta.distributed.cluster.exceptions import (
    trace_protocol_exceptions,
    trace_subprocess_exceptions,
)
from pyrosetta.distributed.cluster.logging_support import (
    bind_protocol,
    setup_target_logging,
    SOCKET_LOGGER_PLUGIN_NAME,
)
from pyrosetta.distributed.cluster.serialization import Serialization
from pyrosetta.distributed.cluster.validators import _validate_residue_type_sets

from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import (
    AbstractSet,
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Tuple,
    TypeVar,
    Union,
)


Q = TypeVar("Q", bound=billiard.Queue)
P = TypeVar("P", bound=billiard.context.Process)
S = TypeVar("S", bound=Serialization)


def _maybe_delay(dt: float, max_delay_time: Union[float, int]) -> None:
    """Maybe delay the user-provided PyRosetta protocol result(s)."""
    delay_time = max_delay_time - dt
    if delay_time > 0.0:
        time.sleep(delay_time)


@trace_protocol_exceptions
def user_protocol(
    packed_pose: PackedPose,
    protocol: Callable[..., Any],
    ignore_errors: bool,
    **kwargs: Dict[Any, Any],
) -> Any:
    """Run the user-provided PyRosetta protocol."""
    with tempfile.TemporaryDirectory() as tmp_path:
        kwargs["PyRosettaCluster_tmp_path"] = tmp_path
        result = protocol(packed_pose, **kwargs)

    return result


@capture_task_metadata
def run_protocol(
    protocol: Callable[..., Any],
    packed_pose: PackedPose,
    datetime_format: str,
    ignore_errors: bool,
    protocols_key: str,
    decoy_ids: List[int],
    serializer: S,
    **kwargs: Dict[Any, Any],
) -> List[Tuple[bytes, bytes]]:
    """Parse the user-provided PyRosetta protocol results."""

    result = user_protocol(packed_pose, protocol, ignore_errors, **kwargs)
    results = _parse_protocol_results(result, kwargs, protocol.__name__, protocols_key, decoy_ids, serializer)

    return results


@trace_subprocess_exceptions
def get_target_results_kwargs(
    q: Q,
    p: P,
    compressed_kwargs: bytes,
    protocol_name: str,
    timeout: Union[float, int],
    ignore_errors: bool,
) -> List[Tuple[Optional[bytes], bytes]]:
    """Get and parse the billiard subprocess results."""

    if p.is_alive():
        return _parse_target_results(
            [obj for obj in iter(q.get(block=True, timeout=timeout))]
        )
    else:
        return [
            (_parse_empty_queue(protocol_name, ignore_errors), compressed_kwargs),
        ]


@bind_protocol
@setup_target_logging
@requires_init
def target(
    protocol: Callable[..., Any],
    compressed_packed_pose: bytes,
    compressed_kwargs: bytes,
    q: Q,
    logging_level: str,
    socket_listener_address: Tuple[str, int],
    datetime_format: str,
    ignore_errors: bool,
    protocols_key: str,
    decoy_ids: List[int],
    compression: Optional[Union[str, bool]],
    client_residue_type_set: AbstractSet[str],
    client_repr: str,
    **pyrosetta_init_kwargs: Dict[str, Any],
) -> None:
    """A wrapper function for a user-provided PyRosetta protocol."""
    serializer = Serialization(compression=compression)
    packed_pose = serializer.decompress_packed_pose(compressed_packed_pose)
    kwargs = serializer.decompress_kwargs(compressed_kwargs)
    kwargs["PyRosettaCluster_client_repr"] = client_repr
    results = run_protocol(
        protocol, packed_pose, datetime_format, ignore_errors, protocols_key, decoy_ids, serializer, **kwargs
    )
    _validate_residue_type_sets(
        _get_residue_type_set(), client_residue_type_set,
    )
    q.put(results)


@bind_protocol
def user_spawn_thread(
    protocol: Callable[..., Any],
    compressed_packed_pose: bytes,
    compressed_kwargs: bytes,
    pyrosetta_init_kwargs: Dict[str, Any],
) -> List[Tuple[Optional[Union[PackedPose, bytes]], Union[Dict[Any, Any], bytes]]]:
    """Generic worker task using the billiard multiprocessing module."""
    t0 = time.time()
    client_repr = repr(get_client())
    worker = get_worker()

    task_args_plugin = worker.plugins[TASK_ARGS_PLUGIN_NAME]
    decoy_ids = task_args_plugin.decoy_ids
    protocols_key = task_args_plugin.protocols_key
    timeout = task_args_plugin.timeout
    ignore_errors = task_args_plugin.ignore_errors
    datetime_format = task_args_plugin.datetime_format
    compression = task_args_plugin.compression
    max_delay_time = task_args_plugin.max_delay_time
    client_residue_type_set = task_args_plugin.client_residue_type_set

    socket_logger_plugin = worker.plugins[SOCKET_LOGGER_PLUGIN_NAME]
    logging_level = socket_logger_plugin.logging_level
    socket_listener_address = (socket_logger_plugin.host, socket_logger_plugin.port)

    q = billiard.Queue()
    p = billiard.context.Process(
        target=target,
        args=(
            protocol,
            compressed_packed_pose,
            compressed_kwargs,
            q,
            logging_level,
            socket_listener_address,
            datetime_format,
            ignore_errors,
            protocols_key,
            decoy_ids,
            compression,
            client_residue_type_set,
            client_repr,
        ),
        kwargs=pyrosetta_init_kwargs,
    )
    p.start()
    results = get_target_results_kwargs(
        q, p, compressed_kwargs, protocol.__name__, timeout, ignore_errors
    )
    p.join()

    dt = time.time() - t0
    _maybe_delay(dt, max_delay_time)

    return results
