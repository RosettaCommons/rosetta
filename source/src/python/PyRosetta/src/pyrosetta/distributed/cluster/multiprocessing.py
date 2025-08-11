# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    import billiard
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.multiprocessing' requires the "
        + "third-party package 'billiard' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/billiard/\n"
    )
    raise

import logging
import tempfile
import time

from pyrosetta.distributed import requires_init
from pyrosetta.distributed.cluster.base import (
    _get_residue_type_set,
    capture_task_metadata,
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
    setup_target_logging,
    setup_worker_logger,
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


def _maybe_delay(dt: float, max_delay_time: Union[float, int], logger: logging.Logger) -> None:
    """Maybe delay the user-provided PyRosetta protocol result(s)."""
    delay_time = max_delay_time - dt
    if delay_time > 0.0:
        logger.info(f"Delaying dask worker results for {delay_time:0.6f} seconds.")
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
    protocol_name: str,
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
    results = _parse_protocol_results(result, kwargs, protocol_name, protocols_key, decoy_ids, serializer)

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
        return _parse_target_results(q.get(block=True, timeout=timeout))
    else:
        return [
            (_parse_empty_queue(protocol_name, ignore_errors), compressed_kwargs),
        ]


@setup_target_logging
@requires_init
def target(
    protocol_name: str,
    compressed_protocol: bytes,
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
    masked_key: bytes,
    task_id: str,
    **pyrosetta_init_kwargs: Dict[str, Any],
) -> None:
    """A wrapper function for a user-provided PyRosetta protocol."""
    del masked_key
    serializer = Serialization(compression=compression)
    protocol = serializer.decompress_object(compressed_protocol)
    packed_pose = serializer.decompress_packed_pose(compressed_packed_pose)
    kwargs = serializer.decompress_kwargs(compressed_kwargs)
    kwargs["PyRosettaCluster_client_repr"] = client_repr
    results = run_protocol(
        protocol_name,
        protocol,
        packed_pose,
        datetime_format,
        ignore_errors,
        protocols_key,
        decoy_ids,
        serializer,
        **kwargs,
    )
    _validate_residue_type_sets(
        _get_residue_type_set(), client_residue_type_set,
    )
    q.put(results)


def user_spawn_thread(
    protocol_name: str,
    compressed_protocol: bytes,
    compressed_packed_pose: bytes,
    compressed_kwargs: bytes,
    pyrosetta_init_kwargs: Dict[str, Any],
    client_repr: str,
    extra_args: Dict[str, Any],
    masked_key: bytes,
    task_id: str,
) -> List[Tuple[Optional[Union[PackedPose, bytes]], Union[Dict[Any, Any], bytes]]]:
    """Generic worker task using the billiard multiprocessing module."""
    t0 = time.time()

    decoy_ids = extra_args["decoy_ids"]
    protocols_key = extra_args["protocols_key"]
    timeout = extra_args["timeout"]
    ignore_errors = extra_args["ignore_errors"]
    datetime_format = extra_args["datetime_format"]
    compression = extra_args["compression"]
    max_delay_time = extra_args["max_delay_time"]
    logging_level = extra_args["logging_level"]
    socket_listener_address = extra_args["socket_listener_address"]
    client_residue_type_set = extra_args["client_residue_type_set"]

    logger = setup_worker_logger(protocol_name, socket_listener_address, masked_key, task_id)
    del masked_key

    # Set the start method to 'spawn' to prevent subprocesses from
    # inheriting PyRosetta's already initialized static singletons
    context = billiard.get_context("spawn")
    q = context.Queue()
    p = context.Process(
        target=target,
        args=(
            protocol_name,
            compressed_protocol,
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
            masked_key,
            task_id,
        ),
        kwargs=pyrosetta_init_kwargs,
    )
    p.start()
    results = get_target_results_kwargs(
        q, p, compressed_kwargs, protocol_name, timeout, ignore_errors
    )
    p.join()

    dt = time.time() - t0
    _maybe_delay(dt, max_delay_time, logger)

    return results
