# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import billiard
    from dask.distributed import get_client
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.multiprocessing' requires the "
        + "third-party packages 'billiard' and 'dask.distributed' as a dependencies!\n"
        + "Please install the package into your python environment. "
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
from pyrosetta.distributed.cluster.logging_support import setup_target_logging
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
    DATETIME_FORMAT: str,
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


@setup_target_logging
@requires_init
def target(
    protocol: Callable[..., Any],
    compressed_packed_pose: bytes,
    compressed_kwargs: bytes,
    q: Q,
    logging_file: str,
    logging_level: str,
    DATETIME_FORMAT: str,
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
        protocol, packed_pose, DATETIME_FORMAT, ignore_errors, protocols_key, decoy_ids, serializer, **kwargs
    )
    _validate_residue_type_sets(
        _get_residue_type_set(), client_residue_type_set,
    )
    q.put(results)


def user_spawn_thread(
    protocol: Callable[..., Any],
    compressed_packed_pose: bytes,
    compressed_kwargs: bytes,
    pyrosetta_init_kwargs: Dict[str, Any],
    decoy_ids: List[int],
    protocols_key: str,
    timeout: Union[float, int],
    ignore_errors: bool,
    logging_file: str,
    logging_level: str,
    DATETIME_FORMAT: str,
    compression: Optional[Union[str, bool]],
    max_delay_time: Union[float, int],
    client_residue_type_set: AbstractSet[str],
) -> List[Tuple[Optional[Union[PackedPose, bytes]], Union[Dict[Any, Any], bytes]]]:
    """Generic worker task using the billiard multiprocessing module."""
    t0 = time.time()
    client_repr = repr(get_client())

    q = billiard.Queue()
    p = billiard.context.Process(
        target=target,
        args=(
            protocol,
            compressed_packed_pose,
            compressed_kwargs,
            q,
            logging_file,
            logging_level,
            DATETIME_FORMAT,
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
