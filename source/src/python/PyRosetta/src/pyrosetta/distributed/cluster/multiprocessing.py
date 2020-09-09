# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"
__email__ = "klima.jason@gmail.com"

try:
    import billiard
    import cloudpickle
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.multiprocessing' requires the "
        + "third-party packages 'billiard' and 'cloudpickle' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/billiard/\n"
        + "https://pypi.org/project/cloudpickle/\n"
    )
    raise

import tempfile

from pyrosetta.distributed import requires_init
from pyrosetta.distributed.cluster.base import (
    _get_packed_pose_kwargs_pairs_list,
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
from pyrosetta.distributed.cluster.validators import _validate_residue_type_sets
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import (
    AbstractSet,
    Any,
    Callable,
    Dict,
    List,
    Tuple,
    TypeVar,
    Union,
)


Q = TypeVar("Q", bound=billiard.Queue)
P = TypeVar("P", bound=billiard.context.Process)


@trace_protocol_exceptions
def user_protocol(
    pose: PackedPose, protocol: Callable[..., Any], ignore_errors: bool, **kwargs: Any
) -> Any:
    """Run the user-provided PyRosetta protocol."""

    with tempfile.TemporaryDirectory():
        result = protocol(pose, **kwargs)

    return result


@capture_task_metadata
def run_protocol(
    protocol: Callable[..., Any],
    pose: PackedPose,
    DATETIME_FORMAT: str,
    ignore_errors: bool,
    **kwargs: Any,
) -> List[Union[PackedPose, bytes]]:
    """Parse the user-provided PyRosetta protocol results."""

    result = user_protocol(pose, protocol, ignore_errors, **kwargs)
    results = _parse_protocol_results(result, protocol.__name__)
    results.append(cloudpickle.dumps(kwargs))

    return results


@trace_subprocess_exceptions
def get_target_results_kwargs(
    q: Q,
    p: P,
    kwargs: Dict[Any, Any],
    protocol_name: str,
    timeout: Union[float, int],
    ignore_errors: bool,
) -> Tuple[List[PackedPose], Dict[Any, Any]]:
    """Get and parse the billiard subprocess results."""

    if p.is_alive():
        return _parse_target_results(
            [obj for obj in iter(q.get(block=True, timeout=timeout))]
        )
    else:
        return _parse_empty_queue(protocol_name), kwargs


@setup_target_logging
@requires_init
def target(
    protocol: Callable[..., Any],
    pose: PackedPose,
    q: Q,
    logging_file: str,
    logging_level: str,
    DATETIME_FORMAT: str,
    ignore_errors: bool,
    master_residue_type_set: AbstractSet[str],
    **kwargs: Any,
) -> None:
    """A wrapper function for a user-provided PyRosetta protocol."""

    results = run_protocol(protocol, pose, DATETIME_FORMAT, ignore_errors, **kwargs,)
    _validate_residue_type_sets(
        _get_residue_type_set(), master_residue_type_set,
    )
    q.put(results)


def user_spawn_thread(
    protocol: Callable[..., Any],
    pose: PackedPose,
    kwargs: Dict[Any, Any],
    decoy_ids: List[int],
    protocols_key: str,
    timeout: Union[float, int],
    ignore_errors: bool,
    logging_file: str,
    logging_level: str,
    DATETIME_FORMAT: str,
    master_residue_type_set: AbstractSet[str],
) -> List[Tuple[PackedPose, Dict[Any, Any]]]:
    """Generic worker task using the billiard multiprocessing module."""

    q = billiard.Queue()
    p = billiard.context.Process(
        target=target,
        args=(
            protocol,
            pose,
            q,
            logging_file,
            logging_level,
            DATETIME_FORMAT,
            ignore_errors,
            master_residue_type_set,
        ),
        kwargs=kwargs,
    )
    p.start()
    results, kwargs = get_target_results_kwargs(
        q, p, kwargs, protocol.__name__, timeout, ignore_errors
    )
    p.join()

    return _get_packed_pose_kwargs_pairs_list(
        results, kwargs, protocol.__name__, protocols_key, decoy_ids
    )
