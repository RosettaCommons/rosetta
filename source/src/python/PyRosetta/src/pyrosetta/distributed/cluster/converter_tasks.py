# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import distributed
    import pandas
    import toolz
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.converter_tasks' requires the "
        + "third-party packages 'distributed', 'pandas', and 'toolz' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/distributed/\n"
        + "https://pypi.org/project/pandas/\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import bz2
import collections
import logging
import json
import os
import pyrosetta.distributed.io as io
import shutil
import subprocess
import warnings

from contextlib import contextmanager
from functools import singledispatch
from pyrosetta.distributed.cluster.config import (
    get_environment_cmd,
    get_environment_manager,
    source_domains,
)
from pyrosetta.distributed.cluster.exceptions import (
    InputError,
    InputFileError,
    OutputError,
)
from pyrosetta.distributed.cluster.io import IO, secure_read_pickle
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.rosetta.basic import was_init_called
from pyrosetta.rosetta.core.pose import Pose
from typing import (
    Any,
    Callable,
    Dict,
    Generator,
    Iterable,
    List,
    NoReturn,
    Optional,
    TypeVar,
    Union,
)


@contextmanager
def not_on_worker() -> Generator[None, Any, None]:
    """A context manager for running code on the host process."""
    try:
        distributed.get_worker()
    except BaseException:
        yield


def maybe_issue_environment_warnings() -> None:
    """
    Issue a warning message if an environment manager is not installed and we are
    not in an active virtual environment on the host process.
    """

    with not_on_worker():
        environment_manager = get_environment_manager()
        if shutil.which(environment_manager):  # An environment manager is installed
            if get_yml() == "":
                warnings.warn(
                    "To use the `pyrosetta.distributed.cluster` namespace and ensure "
                    + "reproducibility of PyRosetta simulations, please either:\n"
                    + "(1) Create and activate a conda or mamba environment (other than 'base'). For instructions, visit:\n"
                    + "https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html\n"
                    + "https://conda.io/activation\n"
                    + "https://mamba.readthedocs.io/en/latest/user_guide/mamba.html\n"
                    + "(2) Create a uv project. For instructions, visit:\n"
                    + "https://docs.astral.sh/uv/getting-started/installation\n"
                    + "https://docs.astral.sh/uv/concepts/projects/init\n"
                    + "(3) Create a pixi manifest. For instructions, visit:\n"
                    + "https://pixi.sh/latest/installation\n"
                    + "https://pixi.sh/latest/getting_started\n",
                    UserWarning,
                    stacklevel=4,
                )  # Warn that we are not in an active virtual environment
        else:  # An environment manager is not installed
            warnings.warn(
                f"The environment manager '{environment_manager}' is not an executable! "
                + "Use of `pyrosetta.distributed.cluster` namespace requires 'conda', 'mamba', "
                + "'uv', or 'pixi' to be properly installed for reproducibility of PyRosetta "
                + "simulations. Please install one of the environment managers onto your system "
                + f"to enable running `which {environment_manager}`. For installation instructions, visit:\n"
                + "https://docs.anaconda.com/anaconda/install\n"
                + "https://github.com/conda-forge/miniforge\n"
                + "https://docs.astral.sh/uv/getting-started/installation\n"
                + "https://pixi.sh/latest/installation\n",
                UserWarning,
                stacklevel=4,
            )  # Warn that environment manager is not in $PATH


def get_protocols_list_of_str(
    input_file: Optional[str] = None,
    scorefile: Optional[str] = None,
    decoy_name: Optional[str] = None,
) -> Union[List[str], NoReturn]:
    """
    Get the user-defined PyRosetta protocols as a `list` object of `str` objects.

    Args:
        input_file: A `str` object specifying the path to the '.pdb' or '.pdb.bz2'
            file from which to extract PyRosettaCluster instance kwargs. If input_file
            is provided, then ignore the scorefile and decoy_name argument parameters.
            Default: None
        scorefile: A `str` object specifying the path to the JSON-formatted scorefile
            from which to extract PyRosettaCluster instance kwargs. If 'scorefile'
            is provided, 'decoy_name' must also be provided. In order to use a scorefile,
            it must contain full simulation records from the original production
            run; i.e., the attribute 'simulation_records_in_scorefile' was set to True.
            Default: None
        decoy_name: A `str` object specifying the decoy name for which to extract
            PyRosettaCluster instance kwargs. If decoy_name is provided, scorefile
            must also be provided.
            Default: None

    Returns:
        A `list` object of `str` objects specifying user-defined PyRosetta protocol names.
    """

    _simulation_records_in_scorefile_msg = (
        "The 'scorefile' parameter argument does not contain the full simulation records. "
        + "In order to reproduce a decoy using a 'scorefile', the PyRosettaCluster "
        + "attribute 'simulation_records_in_scorefile' must have been set to `True` in "
        + "the original simulation. Please provide an 'input_file' generated by PyRosettaCluster, "
        + "or a 'scorefile' with full simulation records generated by PyRosettaCluster, "
        + "in order to reproduce."
    )
    if input_file:
        if scorefile or decoy_name:
            logging.warning(
                "`get_protocols_list_of_str()` received `input_file` and either `scorefile` "
                + " or `decoy_name` keyword argument parameters. Ignoring `scorefile` or "
                + "`decoy_name` keyword argument parameters and using `input_file`!"
            )
        protocols_list_of_str = parse_input_file_to_protocols_str(input_file)
    elif scorefile and decoy_name:
        scorefile = parse_scorefile(scorefile)
        decoy_name = parse_decoy_name(decoy_name)
        if scorefile.endswith(".json"):
            with open(scorefile, "r") as f:
                lines = f.readlines()
                for line in lines:
                    try:
                        scorefile_entry = json.loads(line)
                    except:
                        raise IOError(
                            "`get_protocols_list_of_str()` received `scorefile` which does "
                            + "not appear to be JSON-formatted."
                        )
                    if all(k in scorefile_entry for k in ("metadata", "instance")):
                        if "decoy_name" in scorefile_entry["metadata"]:
                            if scorefile_entry["metadata"]["decoy_name"] == decoy_name:
                                if "protocols" in scorefile_entry["metadata"]:
                                    protocols_list_of_str = scorefile_entry["metadata"]["protocols"]
                                    break
                                else:
                                    raise KeyError(
                                        "'protocols' key not found in 'metadata' entry!"
                                    )
                    else:
                        raise NotImplementedError(_simulation_records_in_scorefile_msg)
        else:
            try:
                df = secure_read_pickle(scorefile, compression="infer")
            except:
                raise TypeError(
                    "`get_protocols_list_of_str()` received `scorefile` which does not appear to be "
                    + "readable by `pyrosetta.distributed.cluster.io.secure_read_pickle(compression='infer')`."
                )
            if all(k in df.columns for k in ("metadata", "instance")):
                for instance, metadata in df[["instance", "metadata"]].values:
                    if "decoy_name" in metadata:
                        if metadata["decoy_name"] == decoy_name:
                            protocols_list_of_str = metadata["protocols"]
                            break
            else:
                raise NotImplementedError(_simulation_records_in_scorefile_msg)
        if not protocols_list_of_str:
            raise KeyError(
                "Error in `get_protocols_list_of_str()`! `decoy_name` is not in `scorefile`."
            )
    else:
        raise NotImplementedError(
            "`get_protocols_list_of_str()` requires either `input_file` or `scorefile` "
            + "and `decoy_name` argument parameter inputs."
        )

    return protocols_list_of_str


def get_scores_dict(obj):
    """Get the PyRosettaCluster scores dictionary from a .pdb or .pdb.bz2 file."""

    if not os.path.exists(obj):
        raise IOError(
            "The `input_file` argument parameter must exist! Received {0}".format(obj)
        )
    else:
        if obj.endswith(".pdb.bz2"):
            with open(obj, "rb") as fbz2:
                pdbstring = bz2.decompress(fbz2.read()).decode()
        elif obj.endswith(".pdb"):
            with open(obj, "r") as f:
                pdbstring = f.read()
        elif obj.endswith((".pkl_pose", ".pkl_pose.bz2", ".b64_pose", ".b64_pose.bz2")):
            if not was_init_called():
                raise RuntimeError(
                    "To get the PyRosettaCluster scores dictionary from a '.pkl_pose', '.pkl_pose.bz2', "
                    + "'.b64_pose' or '.b64_pose.bz2' file, PyRosetta must be initialized (with the same "
                    + "residue type set that was used to save the original file)."
                )
            if obj.endswith(".pkl_pose.bz2"):
                with open(obj, "rb") as fbz2:
                    pdbstring = io.to_pdbstring(io.to_pose(bz2.decompress(fbz2.read())))
            elif obj.endswith(".pkl_pose"):
                with open(obj, "rb") as f:
                    pdbstring = io.to_pdbstring(io.to_pose(f.read()))
            elif obj.endswith(".b64_pose.bz2"):
                with open(obj, "rb") as fbz2:
                    pdbstring = io.to_pdbstring(io.to_pose(bz2.decompress(fbz2.read()).decode()))
            elif obj.endswith(".b64_pose"):
                with open(obj, "r") as f:
                    pdbstring = io.to_pdbstring(io.to_pose(f.read()))
        scores_dict = None
        for line in reversed(pdbstring.split(os.linesep)):
            if line.startswith(IO.REMARK_FORMAT):
                scores_dict = json.loads(line.split(IO.REMARK_FORMAT)[-1])
                break
        else:
            raise IOError(
                "The `input_file` argument parameter must end in "
                + "'.pdb', '.pdb.bz2', '.b64_pose', or '.b64_pose.bz2'."
            )

        if scores_dict is None:
            raise IOError("Could not parse the `input_file` argument parameter!")
        if not all(d in scores_dict for d in ["instance", "metadata", "scores"]):
            raise KeyError("Could not parse the `input_file` argument parameter!")

        return scores_dict


def get_yml() -> str:
    """
    Run environment export command to return a YML file string with the current virtual
    enviroment, excluding certain source domains.
    """

    environment_cmd = get_environment_cmd()
    try:
        raw_yml = subprocess.check_output(
            environment_cmd,
            shell=True,
            stderr=subprocess.DEVNULL,
        ).decode()
    except subprocess.CalledProcessError:
        raw_yml = ""

    return (
        (
            os.linesep.join(
                [
                    line
                    for line in raw_yml.split(os.linesep)
                    if all(
                        source_domain not in line for source_domain in source_domains
                    )
                    and all(not line.startswith(s) for s in ["name:", "prefix:"])
                    and line
                ]
            )
            + os.linesep
        )
        if raw_yml
        else raw_yml
    )


@singledispatch
def to_iterable(obj: Any, func: Callable[..., Any], attr: str) -> List[Any]:
    return [func(obj, attr)]


@to_iterable.register(Pose)
@to_iterable.register(PackedPose)
@to_iterable.register(dict)
def _catch_pose_or_kwargs(
    obj: Union[Pose, PackedPose, Dict[Any, Any]], func: Callable[..., Any], attr: str
) -> List[Any]:
    return [func(obj, attr)]


@to_iterable.register(collections.abc.Iterable)
def _iterate(objs: Iterable[Any], func: Callable[..., Any], attr: str) -> List[Any]:
    return [func(obj, attr) for obj in objs]


@singledispatch
def to_int(obj: Any, attribute: str) -> Union[int, NoReturn]:
    try:
        return int(obj)
    except:
        raise InputError(obj, attribute)


@to_int.register(int)
def _is_int(obj: int, attribute: str) -> int:
    return obj


@singledispatch
def to_packed(obj: Any, protocol_name: str) -> NoReturn:
    """Parse a single result from the user-provided PyRosetta protocol."""

    logging.error(
        f"{protocol_name} did not return objects of type `NoneType`, `Pose`, `PackedPose`, or `dict`!"
    )
    raise OutputError(obj)


@to_packed.register(Pose)
def _to_packed(obj: Pose, protocol_name: str) -> PackedPose:
    return io.to_packed(obj)


@to_packed.register(PackedPose)
@to_packed.register(dict)
def _is_packed_or_kwargs(obj: Union[PackedPose, Dict[Any, Any]], protocol_name: str) -> PackedPose:
    return obj


@to_packed.register(type(None))
def _none_to_packed(obj: None, protocol_name: str) -> PackedPose:
    logging.warning(
        f"{protocol_name} returned `None`. "
        + "Putting an empty `PackedPose` object into the queue."
    )
    return io.to_packed(Pose())


@singledispatch
def to_str(obj: Any, attribute: str) -> Union[str, NoReturn]:
    try:
        return str(int(obj))
    except:
        raise InputError(obj, attribute)


@to_str.register(int)
def _to_int(obj: int, attribute: str) -> str:
    return str(obj)


@to_int.register(float)
@to_str.register(float)
def _to_float(obj: float, attribute: str) -> NoReturn:
    raise NotImplementedError(
        f"PyRosettaCluster '{attribute}' attribute cannot be of type `float`. "
        + f"Received {obj}."
    )


@singledispatch
def parse_input_file_to_protocols_str(obj: Any) -> NoReturn:
    raise InputFileError(obj)


@parse_input_file_to_protocols_str.register(str)
def _parse_str(obj: str) -> List[str]:
    scores_dict = get_scores_dict(obj)
    return scores_dict["metadata"]["protocols"]


@singledispatch
def parse_input_file_to_instance_kwargs(obj: Any) -> NoReturn:
    raise InputFileError(obj)


@parse_input_file_to_instance_kwargs.register(str)
def _parse_str(obj: str) -> Dict[str, Any]:
    scores_dict = get_scores_dict(obj)
    return scores_dict["instance"]


@singledispatch
def parse_scorefile(obj: Any) -> NoReturn:
    raise TypeError(
        "The `scorefile` argument parameter must be of type `str`, "
        + "not of type {0}.".format(type(obj))
    )


@parse_scorefile.register(str)
def _parse_str(obj: str) -> Union[str, NoReturn]:
    if not os.path.exists(obj):
        raise ValueError(
            "The `scorefile` argument parameter must exist! Received {0}".format(obj)
        )
    return obj


@singledispatch
def parse_decoy_name(obj: Any) -> NoReturn:
    raise TypeError(
        "The `decoy_name` argument parameter must be of type `str`, "
        + "not of type {0}.".format(type(obj))
    )


@parse_decoy_name.register(str)
def _from_str(obj: str) -> str:
    return obj


@singledispatch
def reserve_scores_in_results(
    obj: Any, _scores_dict: Dict[Any, Any], protocol_name: str
) -> NoReturn:
    raise OutputError(obj)


@reserve_scores_in_results.register(Pose)
@reserve_scores_in_results.register(PackedPose)
def _parse_packed(
    obj: Union[Pose, PackedPose], _scores_dict: Dict[Any, Any], protocol_name: str
) -> List[PackedPose]:
    packed = to_packed(obj, protocol_name)
    packed.scores = toolz.dicttoolz.merge(_scores_dict, packed.scores)
    return [packed]


@reserve_scores_in_results.register(collections.abc.Iterable)
def _parse_iterable(
    objs: Iterable[Any], _scores_dict: Dict[Any, Any], protocol_name: str
) -> List[PackedPose]:
    out = []
    for obj in objs:
        packed = to_packed(obj, protocol_name)
        if isinstance(packed, PackedPose):
            packed.scores = toolz.dicttoolz.merge(_scores_dict, packed.scores)
        out.append(packed)
    return out


@reserve_scores_in_results.register(type(None))
def _default_none(
    obj: None, _scores_dict: Dict[Any, Any], protocol_name: str
) -> PackedPose:
    return to_packed(obj, protocol_name)


@singledispatch
def parse_client(obj: Any) -> NoReturn:
    raise TypeError(
        "The `client` argument parameter must be of type `distributed.client.Client` "
        + "or `NoneType`, not of type {0}.".format(type(obj))
    )


ClientType = TypeVar("ClientType", bound=distributed.client.Client)


@parse_client.register(distributed.client.Client)
@parse_client.register(type(None))
def _default(obj: Optional[ClientType]) -> Optional[ClientType]:
    return obj


@singledispatch
def parse_input_packed_pose(obj: Any) -> NoReturn:
    raise TypeError(
        "The `input_packed_pose` argument parameter must be of type `PackedPose`, "
        + "`Pose` or `NoneType`, not of type {0}.".format(type(obj))
    )


@parse_input_packed_pose.register(PackedPose)
@parse_input_packed_pose.register(type(None))
def _from_packed_or_none(obj: Optional[PackedPose]) -> Optional[PackedPose]:
    return obj


@parse_input_packed_pose.register(Pose)
def _from_pose(obj: Pose) -> PackedPose:
    return io.to_packed(obj)


@singledispatch
def parse_instance_kwargs(obj: Any) -> NoReturn:
    raise TypeError(
        "The `instance_kwargs` argument parameter must be of type `dict` or "
        + "`NoneType`, not of type {0}.".format(type(obj))
    )


@parse_instance_kwargs.register(dict)
def _parse_dict(obj: Dict[Any, Any]) -> Dict[Any, Any]:
    for k in obj.keys():
        if k in ("client", "clients", "input_packed_pose"):
            raise NotImplementedError(
                f"The parameter '{k}' must be passed directly to `reproduce()`, "
                + "not as a member of the 'instance_kwargs' dictionary."
            )
        elif k in ("seeds", "decoy_ids"):
            raise NotImplementedError(
                f"The parameter '{k}' must be obtained from the original input file "
                + "or scorefile, not input as a member of the 'instance_kwargs' dictionary."
            )
        elif k == "filter_results":
            raise ValueError(
                f"The parameter '{k}' cannot be set as a PyRosettaCluster attribute "
                + "in `reproduce()` because the saved 'decoy_ids' attribute from the "
                + "original simulation depends on the original decoy output order from "
                + "each protocol, so results must be filtered identically. Please remove "
                + "this keyword argument to run `reproduce()`."
            )
    return obj


@parse_instance_kwargs.register(type(None))
def _default_none(obj: None) -> Dict[Any, Any]:
    return {}


@singledispatch
def is_empty(obj: Any) -> NoReturn:
    """Test whether a `PackedPose` object is empty."""
    raise NotImplementedError(type(obj))

@is_empty.register(type(None))
def _from_none(obj: None) -> bool:
    # Protocol results return a `None` object when a segmentation fault occurs with `ignore_errors=True`
    return False

@is_empty.register(PackedPose)
def _from_packed(obj: PackedPose) -> bool:
    return obj.empty()


def is_bytes(obj: Any) -> bool:
    return isinstance(obj, bytes)


def is_packed(obj: Any) -> bool:
    return isinstance(obj, PackedPose)


def is_dict(obj: Any) -> bool:
    return isinstance(obj, dict)
