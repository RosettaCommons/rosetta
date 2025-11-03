# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import distributed
    import toolz
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.tools' requires the "
        + "third-party packages 'distributed' and 'toolz' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/distributed/\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import inspect
import json
import logging
import os
import subprocess
import sys
import tempfile
import warnings

from datetime import datetime
from functools import wraps
from pyrosetta.distributed.cluster.config import get_environment_config
from pyrosetta.distributed.cluster.converters import _parse_protocols
from pyrosetta.distributed.cluster.converter_tasks import (
    is_dict,
    is_empty,
    get_protocols_list_of_str,
    parse_client,
    parse_decoy_name,
    parse_init_file,
    parse_input_file_to_instance_kwargs,
    parse_instance_kwargs,
    parse_scorefile,
    reserve_scores_in_results,
)
from pyrosetta.distributed.cluster.io import secure_read_pickle
from pyrosetta.distributed.cluster.serialization import (
    Serialization,
    update_scores,
)
from pyrosetta.distributed.cluster.core import PyRosettaCluster
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.rosetta.basic import was_init_called
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.exceptions import PyRosettaIsNotInitializedError
from typing import (
    Any,
    Callable,
    Dict,
    Generator,
    List,
    NoReturn,
    Optional,
    Tuple,
    TypeVar,
    Union,
    cast,
)


P = TypeVar("P", bound=Callable[..., Any])


def get_protocols(
    protocols: Optional[Union[List[Union[Callable[..., Any], str]], Callable[..., Any], str]] = None,
    input_file: Optional[Union[str, Pose, PackedPose]] = None,
    scorefile: Optional[str] = None,
    decoy_name: Optional[str] = None,
) -> Union[List[Union[Callable[..., Any], str]], NoReturn]:
    """
    Given an 'input_file' that was written by PyRosettaCluster, or a full 'scorefile'
    and a 'decoy_name' that was written by PyRosettaCluster, if 'protocols' is provided
    then validate the 'protocols' against those in the 'input_file' or 'scorefile',
    otherwise if 'protocols' is `NoneType` then attempt to return the PyRosettaCluster
    protocols from the current scope matching the protocol names in the 'input_file'
    or 'scorefile' keyword argument parameters.

    Args:
        protocols: An iterable of `str` objects specifying the names of user-provided
            PyRosetta protocols to validate or return.
            Default: None
        input_file: A `str` object specifying the path to the '.pdb', '.pdb.bz2', '.pkl_pose',
            '.pkl_pose.bz2', '.b64_pose', '.b64_pose.bz2', '.init', or '.init.bz2' file,
            or a `Pose`or `PackedPose` object, from which to extract PyRosettaCluster instance
            kwargs. If 'input_file' is provided, then ignore the 'scorefile' and 'decoy_name'
            keyword argument parameters.
            Default: None
        scorefile: A `str` object specifying the path to the JSON-formatted scorefile
            (or pickled `pandas.DataFrame` scorefile) from a PyRosettaCluster simulation
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
        A `list` of user-defined PyRosetta protocol names from the 'input_file' or 'scorefile'.
        If `protocols` is None, then attempt to return the PyRosettaCluster protocols
        from the current scope matching the protocol names in the 'input_file' or 'scorefile'.
    """

    if protocols:
        # Validate the user-provided protocols against the original list of protocol name strings
        input_protocols_list_of_str = [
            protocol.__name__ for protocol in _parse_protocols(protocols)
        ]
        original_protocols_list_of_str = get_protocols_list_of_str(
            input_file=input_file, scorefile=scorefile, decoy_name=decoy_name
        )
        assert len(original_protocols_list_of_str) == len(
            input_protocols_list_of_str
        ), (
            "The original user-defined PyRosetta protocols list and the 'protocols' argument "
            + " parameter have different lengths! Cannot reproduce!"
        )
        if original_protocols_list_of_str == input_protocols_list_of_str:
            logging.info(
                "The 'protocols' argument parameter matches the user-defined PyRosetta protocols "
                + "from the original production run. Continuing with the reproduction."
            )
        for i in range(len(original_protocols_list_of_str)):
            if original_protocols_list_of_str[i] != input_protocols_list_of_str[i]:
                logging.warning(
                    f"The original user-defined PyRosetta protocol '{original_protocols_list_of_str[i]}' "
                    + f"appears to have changed names to '{input_protocols_list_of_str[i]}'! "
                    + f"Please verify that {input_protocols_list_of_str[i]}' is functionally equivalent to "
                    + f"'{original_protocols_list_of_str[i]}', otherwise the reproduction simulation "
                    + "will not reproduce the original decoy(s)! Continuing with the reproduction run."
                )
    else:
        # Get the protocols in the scope from the list of protocol names as strings
        # Use locals and globals back two frames to the scope of the caller of `reproduce`
        scope = inspect.currentframe().f_back.f_back
        protocols = []
        for protocol_name in get_protocols_list_of_str(
            input_file=input_file, scorefile=scorefile, decoy_name=decoy_name
        ):
            if protocol_name in scope.f_locals:
                logging.info(
                    f"Automatically detected user-provided PyRosetta protocol '{protocol_name}' "
                    + "in the local variables of the current frame."
                )
                protocols.append(scope.f_locals[protocol_name])
            elif protocol_name in scope.f_globals:
                logging.info(
                    f"Automatically detected user-provided PyRosetta protocol '{protocol_name}' "
                    + "in the global variables of the current frame."
                )
                protocols.append(scope.f_globals[protocol_name])
            else:
                raise RuntimeError(
                    f"The original user-defined PyRosetta protocol '{protocol_name}' "
                    + "could not be found in the current scope! Please verify that the original "
                    + "user-defined PyRosetta protocols are defined in the same scope as they were in "
                    + "the original production run. Alternatively, you may pass the original "
                    + "user-defined PyRosetta protocols (defined under the new scope) maintaining "
                    + "the original protocol order as a parameter to the 'protocols' keyword argument."
                )

    return _parse_protocols(protocols)


def get_instance_kwargs(
    input_file: Optional[Union[str, Pose, PackedPose]] = None,
    scorefile: Optional[str] = None,
    decoy_name: Optional[str] = None,
    skip_corrections: Optional[bool] = None,
) -> Union[Dict[str, Any], NoReturn]:
    """
    Given an input file that was written by PyRosettaCluster, or a scorefile
    and a decoy name that was written by PyRosettaCluster, return the PyRosettaCluster
    instance kwargs needed to reproduce the decoy using PyRosettaCluster.

    Args:
        input_file: A `str` object specifying the path to the '.pdb', '.pdb.bz2', '.pkl_pose',
            '.pkl_pose.bz2', '.b64_pose', '.b64_pose.bz2', '.init', or '.init.bz2' file, or a
            `Pose`or `PackedPose` object, from which to extract PyRosettaCluster instance kwargs.
            If 'input_file' is provided, then ignore the 'scorefile' and 'decoy_name' keyword
            argument parameters.
            Default: None
        scorefile: A `str` object specifying the path to the JSON-formatted scorefile
            (or pickled `pandas.DataFrame` scorefile) from a PyRosettaCluster simulation
            from which to extract PyRosettaCluster instance kwargs. If 'scorefile'
            is provided, 'decoy_name' must also be provided. In order to use a scorefile,
            it must contain full simulation records from the original production
            run; i.e., the attribute 'simulation_records_in_scorefile' was set to True.
            Default: None
        decoy_name: A `str` object specifying the decoy name for which to extract
            PyRosettaCluster instance kwargs. If 'decoy_name' is provided, 'scorefile'
            must also be provided.
            Default: None
        skip_corrections: A `bool` object specifying whether or not to skip any ScoreFunction
            corrections specified in the PyRosettaCluster task initialization options
            (extracted from the 'input_file' or 'scorefile' keyword argument parameter).
            Default: None

    Returns:
        A `dict` object of PyRosettaCluster instance kwargs.
    """
    _simulation_records_in_scorefile_msg = (
        "The 'scorefile' argument parameter does not contain the full simulation records. "
        + "In order to reproduce a decoy using a 'scorefile', the PyRosettaCluster "
        + "attribute 'simulation_records_in_scorefile' must have been set to `True` in "
        + "the original simulation. Please provide an 'input_file' generated by PyRosettaCluster, "
        + "or a 'scorefile' with full simulation records generated by PyRosettaCluster, "
        + "in order to reproduce."
    )
    if input_file:
        if scorefile or decoy_name:
            warnings.warn(
                "Received 'input_file' and either 'scorefile' or 'decoy_name' keyword argument parameters. "
                + "Ignoring 'scorefile' and 'decoy_name' and using 'input_file' keyword argument parameter!",
                UserWarning,
                stacklevel=2,
            )
        instance_kwargs = parse_input_file_to_instance_kwargs(input_file)
    elif scorefile and decoy_name:
        scorefile = parse_scorefile(scorefile)
        decoy_name = parse_decoy_name(decoy_name)
        instance_kwargs = None
        if scorefile.endswith(".json"):
            with open(scorefile, "r") as f:
                lines = f.readlines()
                for line in lines:
                    try:
                        scorefile_entry = json.loads(line)
                    except:
                        raise TypeError(
                            "`get_instance_kwargs()` received `scorefile` which does not appear to be JSON-formatted."
                        )
                    if all(k in scorefile_entry for k in ("metadata", "instance")):
                        if "decoy_name" in scorefile_entry["metadata"]:
                            if scorefile_entry["metadata"]["decoy_name"] == decoy_name:
                                instance_kwargs = scorefile_entry["instance"]
                                break
                    else:
                        raise NotImplementedError(_simulation_records_in_scorefile_msg)
        else:
            try:
                df = secure_read_pickle(scorefile, compression="infer")
            except:
                raise TypeError(
                    "`get_instance_kwargs()` received `scorefile` which does not appear to be "
                    + "readable by `pyrosetta.distributed.cluster.io.secure_read_pickle(compression='infer')`."
                )
            if all(k in df.columns for k in ("metadata", "instance")):
                for instance, metadata in df[["instance", "metadata"]].values:
                    if "decoy_name" in metadata:
                        if metadata["decoy_name"] == decoy_name:
                            instance_kwargs = dict(instance)
                            break
            else:
                raise NotImplementedError(_simulation_records_in_scorefile_msg)
        if instance_kwargs is None:
            raise KeyError(
                "Error in `get_instance_kwargs()`! The provided `decoy_name` is not in the provided `scorefile`."
            )
    else:
        raise NotImplementedError(
            "`get_instance_kwargs()` requires either `input_file` (or `scorefile` and `decoy_name`) argument parameter inputs."
        )
    assert isinstance(
        instance_kwargs, dict
    ), "Returned instance keyword arguments are not of type `dict`."

    if skip_corrections:
        assert isinstance(
            instance_kwargs["tasks"], dict
        ), "PyRosettaCluster 'tasks' keyword argument parameter must be an instance of `dict`."
        for option in ("extra_options", "options"):
            if option in instance_kwargs["tasks"]:
                if isinstance(instance_kwargs["tasks"][option], dict):
                    instance_kwargs["tasks"][option] = toolz.dicttoolz.keyfilter(
                        lambda k: not k.startswith(("corrections:", "-corrections:")),
                        instance_kwargs["tasks"][option],
                    )
                elif isinstance(instance_kwargs["tasks"][option], str):
                    if "corrections:" in instance_kwargs["tasks"][option]:
                        raise NotImplementedError(
                            "Cannot skip ScoreFunction corrections because the original simulation did not output "
                            + "PyRosettaCluster results with normalized PyRosetta initialization options or configure "
                            + "the task's PyRosetta initialization options as an instance of `dict`. Please disable the "
                            + "'skip_corrections' keyword argument to continue with the reproduction."
                        )
                else:
                    raise TypeError(
                        f"PyRosettaCluster task key '{option}' must have a value of type `dict` or `str`. "
                        + f"Received: {type(instance_kwargs['tasks'][option])}"
                    )

    return instance_kwargs


def recreate_environment(
    environment_name: Optional[str] = None,
    input_file: Optional[Union[str, Pose, PackedPose]] = None,
    scorefile: Optional[str] = None,
    decoy_name: Optional[str] = None,
    timeout: Optional[int] = None,
) -> Optional[NoReturn]:
    """
    Given an input file that was written by PyRosettaCluster, or a scorefile
    and a decoy name that was written by PyRosettaCluster, recreate the
    environment that was used to generate the decoy with a new environment name.

    The environment manager used (i.e., either 'conda', 'mamba', 'uv', or 'pixi') is
    automatically determined from the operating system environment variable
    'PYROSETTACLUSTER_ENVIRONMENT_MANAGER' if exported, or otherwise it is automatically
    selected based on which manager is available on the system. If no executables
    are explicitly found, then 'conda' is used by default.

    Args:
        environment_name: A `str` object specifying the new name of the environment
            to recreate. If using 'conda' and 'mamba', this is the prefix directory that will
            be created in the current working directory. If using 'uv' or 'pixi', this is the
            local project directory name that will be created in the current working directory.
            Default: 'PyRosettaCluster_' + datetime.now().strftime("%Y.%m.%d.%H.%M.%S.%f")
        input_file: A `str` object specifying the path to the '.pdb', '.pdb.bz2', '.pkl_pose',
            '.pkl_pose.bz2', '.b64_pose', or '.b64_pose.bz2' file, or a `Pose` or `PackedPose`
            object, from which to extract PyRosettaCluster instance kwargs. If 'input_file' is
            provided, then ignore the 'scorefile' and 'decoy_name' keyword argument parameters.
            Default: None
        scorefile: A `str` object specifying the path to the JSON-formatted scorefile
            (or pickled `pandas.DataFrame` scorefile) from a PyRosettaCluster simulation
            from which to extract PyRosettaCluster instance kwargs. If 'scorefile'
            is provided, 'decoy_name' must also be provided. The scorefile must
            contain full simulation records from the original production run
            (i.e., 'simulation_records_in_scorefile' was set to `True`).
            Default: None
        decoy_name: A `str` object specifying the decoy name for which to extract
            PyRosettaCluster instance kwargs. Must be provided if 'scorefile'
            is used.
            Default: None
        timeout: An `int` object specifying the timeout in seconds before any
            subprocesses are terminated.
            Default: None

    Returns:
        None
    """
    def _run_subprocess(cmd) -> str:
        try:
            return subprocess.check_output(
                cmd,
                shell=True,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                text=True,
            )
        except subprocess.CalledProcessError as ex:
            raise RuntimeError(
                f"Command failed: `{cmd}`\n"
                f"Return code: {ex.returncode}\n"
                f"Output:\n{ex.output}"
            ) from ex

    if not environment_name:
        environment_name = "PyRosettaCluster_" + datetime.now().strftime(
            "%Y.%m.%d.%H.%M.%S.%f"
        )

    _env_config = get_environment_config()
    environment_manager = _env_config.environment_manager
    env_list_cmd = _env_config.env_list_cmd

    # Test if environment name already exists
    if environment_manager in ("conda", "mamba"):
        env_list_cmd = _env_config.env_list_cmd
        envs_out = _run_subprocess(env_list_cmd)
        envs_all = [
            line.split()[0]
            for line in envs_out.splitlines()
            if line and not line.startswith("#")
        ]
    elif environment_manager in ("uv", "pixi"):
        # Only consider local uv or pixi projects created in the current working directory
        if environment_manager == "uv":
            _filenames = ("uv.lock", "pyproject.toml",)
        elif environment_manager == "pixi":
            _filenames = ("pixi.toml",)
        cwd = os.getcwd()
        envs_all = [
            d for d in os.listdir(cwd)
            if os.path.isdir(d) and any(os.path.isfile(os.path.join(cwd, d, f)) for f in _filenames)
        ]
    else:
        raise RuntimeError(f"Unsupported environment manager: '{environment_manager}'")
    if environment_name in envs_all:
        raise RuntimeError(
            f"The environment name '{environment_name}' already exists under {environment_manager}."
        )

    # Extract environment spec from record
    _instance_kwargs = get_instance_kwargs(
        input_file=input_file,
        scorefile=scorefile,
        decoy_name=decoy_name,
        skip_corrections=None,
    )
    if "environment" not in _instance_kwargs:
        raise RuntimeError(
            "PyRosettaCluster 'environment' instance attribute does not exist. "
            + "`recreate_environment()` cannot create environment!"
        )
    raw_spec = _instance_kwargs.get("environment", None)
    if not raw_spec:
        raise RuntimeError(
            "PyRosettaCluster 'environment' instance attribute is empty. "
            + "`recreate_environment()` cannot create environment!"
        )
    # Get original environment manager
    for line in raw_spec.splitlines():
        if line.startswith(f"# {_env_config._ENV_VAR}"):
            original_environment_manager = line.split("=")[-1].strip()
            break
    else:
        original_environment_manager = "conda" # For legacy PyRosettaCluster simulations with a non-empty raw spec

    # Test that current environment manager can recreate original environment
    if original_environment_manager == "uv" and environment_manager != "uv":
        raise RuntimeError(
            f"The original PyRosettaCluster simulation used '{original_environment_manager}' as "
            + f"an environment manager, but the current environment manager is configured to use "
            + f"'{environment_manager}' as an environment manager! The environment specification "
            + f"is saved in a 'requirements.txt' format, and therefore requires '{original_environment_manager}' "
            + f"to recreate the environment. Please ensure '{original_environment_manager}' is installed "
            + f"and run `export {_env_config._ENV_VAR}={original_environment_manager}` to properly configure "
            + f"the '{original_environment_manager}' environment manager, then try again. For installation "
            + "instructions, please visit:\n"
            + "https://docs.astral.sh/uv/guides/install-python\n"
        )
    elif original_environment_manager != "uv" and environment_manager == "uv":
        raise RuntimeError(
            f"The original PyRosettaCluster simulation used '{original_environment_manager}' as "
            + "an environment manager, but the current environment manager is configured to use "
            + f"'{environment_manager}' as an environment manager! The environment specification "
            + "is saved in a YAML file format, and therefore requires 'conda', 'mamba', or 'pixi' to "
            + "recreate the environment. Please ensure that 'conda', 'mamba', or 'pixi' is installed, then "
            + "configure the environment manager by running one of the following commands, then try again:\n"
            + f"To configure 'conda', run: `export {_env_config._ENV_VAR}=conda`\n"
            + f"To configure 'mamba', run: `export {_env_config._ENV_VAR}=mamba`\n"
            + f"To configure 'pixi', run:  `export {_env_config._ENV_VAR}=pixi`\n"
            + "For installation instructions, please visit:\n"
            + "https://docs.anaconda.com/anaconda/install\n"
            + "https://github.com/conda-forge/miniforge\n"
            + "https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html\n"
            + "https://pixi.sh/latest/installation\n"
        )
    elif original_environment_manager != environment_manager:
        warnings.warn(
            f"The original PyRosettaCluster simulation used '{original_environment_manager}' as "
            + "an environment manager, but the current environment manager is configured to use "
            + f"'{environment_manager}' as an environment manager! Now attempting to recreate the "
            + "environment with cross-manager compatibility using the universal YAML file string...",
            RuntimeWarning,
            stacklevel=2,
        )

    # Recreate the environment
    with tempfile.TemporaryDirectory() as tmp_dir:
        env_create_cmd = _env_config.env_create_cmd(environment_name, raw_spec, tmp_dir)
        output = _run_subprocess(env_create_cmd)
        sys.stdout.write(
            f"\nEnvironment successfully created using {environment_manager}: '{environment_name}'\n"
            f"Output:\n{output}\n"
        )


def reserve_scores(func: P) -> Union[P, NoReturn]:
    """
    Use this as a Python decorator of any user-provided PyRosetta protocol.
    If any scoreterms and values are present in the input `packed_pose`, then if
    they are deleted during execution of the decorated user-provided PyRosetta
    protocol, then append those scoreterms and values back into the `pose.cache`
    dictionary after execution. If any scoreterms and values are present in the
    input `packed_pose` and also present in the returned or yielded output `Pose`
    or `PackedPose` objects, then do not append the original scoreterms and values
    back into the `pose.cache` dictionary after execution (that is, keep the outputted
    scoreterms and values in the `pose.cache` dictionary). Any new scoreterms and
    values acquired in the decorated user-provided PyRosetta protocol will never
    be overwritten. This allows users to maintain scoreterms and values acquired
    in earlier user-defined PyRosetta protocols if needing to execute Rosetta
    Movers that happen to delete scores from pose objects.

    For example:

    @reserve_scores
    def my_pyrosetta_protocol(packed_pose, **kwargs):
        from pyrosetta import MyMover
        pose = packed_pose.pose
        MyMover().apply(pose)
        return pose

    Args:
        A user-provided PyRosetta function.

    Returns:
        The output from the user-provided PyRosetta function, reserving the scores.
    """
    import pyrosetta  # noqa
    import pyrosetta.distributed  # noqa

    @wraps(func)
    def wrapper(packed_pose, **kwargs):
        if packed_pose is not None:
            _scores_dict = update_scores(packed_pose).scores
        else:
            _scores_dict = {}
        _output = func(packed_pose, **kwargs)

        return reserve_scores_in_results(_output, _scores_dict, func.__name__)

    return cast(P, wrapper)


def requires_packed_pose(func: P) -> Union[PackedPose, None, P]:
    """
    Use this as a Python decorator of any user-provided PyRosetta protocol.
    If a user-provided PyRosetta protocol requires that the first argument
    parameter be a non-empty `PackedPose` object, then return any received empty
    `PackedPose` objects or `NoneType` objects and skip the decorated protocol,
    otherwise run the decorated protocol.

    If using `PyRosettaCluster(filter_results=False)` and the preceding protocol
    returns or yields either `None`, an empty `Pose` object, or an empty `PackedPose`
    object, then an empty `PackedPose` object is distributed to the next user-provided
    PyRosetta protocol, in which case the next protocol and/or any downstream
    protocols are skipped if they are decorated with this decorator. If using
    `PyRosettaCluster(ignore_errors=True)` and an error is raised in the preceding
    protocol, then a `NoneType` object is distributed to the next user-provided
    PyRosetta protocol, in which case the next protocol and/or any downstream
    protocols are skipped if they are decorated with this decorator.

    For example:

    @requires_packed_pose
    def my_pyrosetta_protocol(packed_pose, **kwargs):
        assert packed_pose.pose.size() > 0
        return packed_pose

    Args:
        A user-provided PyRosetta function.

    Returns:
        The input `packed_pose` argument parameter if it is an empty `PackedPose` object
        or a `NoneType` object, otherwise the results from the decorated protocol.
    """
    @wraps(func)
    def wrapper(packed_pose, **kwargs):
        _msg = "User-provided PyRosetta protocol '{0}' received and is duly returning {1} object."
        if is_empty(packed_pose):
            logging.info(_msg.format(func.__name__, "an empty `PackedPose`"))
            return packed_pose
        elif packed_pose is None:
            logging.info(_msg.format(func.__name__, "a `NoneType`"))
            return packed_pose
        else:
            return func(packed_pose, **kwargs)

    return cast(P, wrapper)


def reproduce(
    input_file: Optional[str] = None,
    scorefile: Optional[str] = None,
    decoy_name: Optional[str] = None,
    protocols: Any = None,
    client: Optional[distributed.client.Client] = None,
    clients: Optional[List[distributed.client.Client]] = None,
    input_packed_pose: Optional[Union[Pose, PackedPose]] = None,
    instance_kwargs: Optional[Dict[Any, Any]] = None,
    clients_indices: Optional[List[int]] = None,
    resources: Optional[Dict[Any, Any]] = None,
    skip_corrections: bool = False,
    init_from_file_kwargs: Optional[Dict[str, Any]] = None,
) -> Optional[NoReturn]:
    """
    Given an input file that was written by PyRosettaCluster (or a full scorefile
    and a decoy name that was written by PyRosettaCluster) and any additional
    PyRosettaCluster instance kwargs, run the reproduction simulation for the
    given decoy with a new instance of PyRosettaCluster.

    Args:
        input_file: A `str` object specifying the path to the '.pdb', '.pdb.bz2',
            '.pkl_pose', '.pkl_pose.bz2', '.b64_pose', '.b64_pose.bz2', '.init' or '.init.bz2'
            file from which to extract PyRosettaCluster instance kwargs. If 'input_file'
            is provided, then ignore the 'scorefile' and 'decoy_name' argument parameters.
            If a '.init' or '.init.bz2' file is provided and PyRosetta is not yet initialized,
            this first initializes PyRosetta with the PyRosetta initialization file (see the
            'init_from_file_kwargs' keyword argument). Note that '.pkl_pose', '.pkl_pose.bz2',
            '.b64_pose', '.b64_pose.bz2', '.init' and '.init.bz2' files contain pickled Pose
            objects that are deserialized using PyRosetta's secure unpickler upon running the
            `reproduce()` function, but please still only input these file types if you know and
            trust their source. Learn more `here <https://docs.python.org/3/library/pickle.html>`_.
            Default: None
        scorefile: A `str` object specifying the path to the JSON-formatted scorefile
            (or pickled `pandas.DataFrame` scorefile) from a PyRosettaCluster simulation
            from which to extract PyRosettaCluster instance kwargs. If 'scorefile'
            is provided, 'decoy_name' must also be provided. In order to use a scorefile,
            it must contain full simulation records from the original production
            run; i.e., the attribute 'simulation_records_in_scorefile' was set to True.
            Note that in order to securely load pickled `pandas.DataFrame` objects, please
            ensure that `pyrosetta.secure_unpickle.add_secure_package("pandas")` has been run.
            Default: None
        decoy_name: A `str` object specifying the decoy name for which to extract
            PyRosettaCluster instance kwargs. If decoy_name is provided, scorefile
            must also be provided.
            Default: None
        protocols: An optional iterable object of function or generator objects specifying
            an ordered sequence of user-defined PyRosetta protocols to execute for
            the reproduction. This argument only needs to be provided if the user-defined
            PyRosetta protocols are not defined with the same scope as in the original
            production run.
            Default: None
        client: An optional initialized dask `distributed.client.Client` object to be used as
            the dask client interface to the local or remote compute cluster. If `None`,
            then PyRosettaCluster initializes its own dask client based on the settings
            from the original production run. Deprecated by the `clients` attribute, but
            supported for legacy purposes.
            Default: None
        clients: A `list` or `tuple` object of initialized dask `distributed.client.Client`
            objects to be used as the dask client interface(s) to the local or remote compute
            cluster(s). If `None`, then PyRosettaCluster initializes its own dask client based
            on the settings from the original production run. Optionally used in
            combination with the `clients_indices` attribute.
            Default: None
        input_packed_pose: An optional input `PackedPose` object that is accessible via
            the first argument of the first user-defined PyRosetta protocol.
            Default: None
        instance_kwargs: An optional `dict` object of valid PyRosettaCluster attributes
            which will override any PyRosettaCluster attributes that were used to generate
            the original decoy.
            Default: None
        clients_indices: An optional `list` or `tuple` object of `int` objects, where each `int` object represents
            a zero-based index corresponding to the initialized dask `distributed.client.Client` object(s) passed 
            to the `PyRosettaCluster(clients=...)` class attribute. If not `None`, then the length of the 
            `clients_indices` object must equal the number of protocols passed to the `PyRosettaCluster().distribute`
            method.
            Default: None
        resources: An optional `list` or `tuple` object of `dict` objects, where each `dict` object represents
            an abstract, arbitrary resource to constrain which dask workers run the user-defined PyRosetta protocols.
            If `None`, then do not impose resource constaints on any protocols. If not `None`, then the length
            of the `resources` object must equal the number of protocols passed to the `PyRosettaCluster().distribute`
            method, such that each resource specified indicates the unique resource constraints for the protocol at the
            corresponding index of the protocols passed to `PyRosettaCluster().distribute`. Note that this feature is only 
            useful when one passes in their own instantiated client(s) with dask workers set up with various resource
            constraints. If dask workers were not instantiated to satisfy the specified resource constraints, protocols
            will hang indefinitely because the dask scheduler is waiting for workers that meet the specified resource 
            constraints so that it can schedule these protocols. Unless workers were created with these resource tags
            applied, the protocols will not run. See https://distributed.dask.org/en/latest/resources.html for more
            information.
            Default: None
        skip_corrections: A `bool` object specifying whether or not to skip any ScoreFunction corrections specified in
            the PyRosettaCluster task 'options' or 'extra_options' values (extracted from either the 'input_file' or
            'scorefile' keyword argument parameter), which are set in-code upon PyRosetta initialization. If the current
            PyRosetta build and conda environment are identical to those used for the original simulation, this parameter
            may be set to `True` to enable the reproduced decoy output file to be used for successive reproductions. If
            reproducing from a '.init' file, it is recommended to also set 'skip_corrections' of the 'init_from_file_kwargs'
            keyword argument to the same value.
            Default: False
        init_from_file_kwargs: An optional `dict` object to override the default `pyrosetta.init_from_file()` keyword
            arguments if the 'input_file' keyword argument parameter is a path to a '.init' file, otherwise it is not used.
            See the `pyrosetta.init_from_file` docstring for more information.
            Default: {
                'output_dir': os.path.join(tempfile.TemporaryDirectory().name, "pyrosetta_init_input_files"),
                'skip_corrections': skip_corrections, # Defaults to the 'skip_corrections' value from `reproduce()`
                'relative_paths': True,
                'dry_run': False,
                'max_decompressed_bytes': pow(2, 30), # 1 GiB
                'database': None,
                'verbose': True,
                'set_logging_handler': 'logging',
                'notebook': None,
                'silent': False,
            }

    Returns:
        None
    """
    if not isinstance(skip_corrections, bool):
        raise TypeError(
            "The 'skip_corrections' keyword argument parameter must be of type `bool`. "
            + f"Received: {type(skip_corrections)}"
        )

    _tmp_dir = None
    if isinstance(input_file, str):
        if input_file.endswith((".init", ".init.bz2")):
            _tmp_dir = tempfile.TemporaryDirectory(prefix="PyRosettaCluster_reproduce_")
            default_init_from_file_kwargs = dict(
                output_dir=os.path.join(_tmp_dir.name, "pyrosetta_init_input_files"),
                skip_corrections=skip_corrections,
                relative_paths=True,
                dry_run=False,
                max_decompressed_bytes=pow(2, 30), # 1 GiB
                database=None,
                verbose=True,
                set_logging_handler="logging",
                notebook=None,
                silent=False,
            )
            input_packed_pose, input_file = parse_init_file(
                input_file,
                input_packed_pose,
                skip_corrections,
                toolz.dicttoolz.merge(
                    default_init_from_file_kwargs,
                    init_from_file_kwargs if is_dict(init_from_file_kwargs) else {},
                ),
            )
        elif input_file.endswith((".pkl_pose", ".pkl_pose.bz2", ".b64_pose", ".b64_pose.bz2")):
            if not was_init_called():
                raise PyRosettaIsNotInitializedError(
                    "If providing a '.pkl_pose', '.pkl_pose.bz2', '.b64_pose', or '.b64_pose.bz2' file to the 'input_file' "
                    + "keyword argument parameter, please ensure `pyrosetta.init()` or `pyrosetta.init_from_file()` has been "
                    + "properly called (with the same residue type set as that used to generate the original '.pkl_pose', "
                    + "'.pkl_pose.bz2', '.b64_pose', or '.b64_pose.bz2' file) before running the `reproduce()` function. "
                    + "If an output '.init' file from the original simulation is available, it is recommended to run "
                    + "`pyrosetta.init_from_file()` with that '.init' file before running the `reproduce()` function."
                )

    PyRosettaCluster(
        **toolz.dicttoolz.keyfilter(
            lambda a: a not in ["client", "clients", "input_packed_pose"],
            toolz.dicttoolz.merge(
                get_instance_kwargs(
                    input_file=input_file,
                    scorefile=scorefile,
                    decoy_name=decoy_name,
                    skip_corrections=skip_corrections,
                ),
                parse_instance_kwargs(instance_kwargs),
            ),
        ),
        client=parse_client(client),
        clients=clients,
        input_packed_pose=input_packed_pose,
    ).distribute(
        protocols=get_protocols(
            protocols=protocols,
            input_file=input_file,
            scorefile=scorefile,
            decoy_name=decoy_name,
        ),
        clients_indices=clients_indices,
        resources=resources,
    )
    if isinstance(_tmp_dir, tempfile.TemporaryDirectory):
        _tmp_dir.cleanup()


def produce(**kwargs: Any) -> Optional[NoReturn]:
    """
    `PyRosettaCluster().distribute()` shim requiring the 'protocols' keyword argument, and optionally
    any PyRosettaCluster keyword arguments or the 'clients_indices' keyword argument (when using
    the `PyRosettaCluster(clients=...)` keyword argument), or the 'resources' keyword argument.

    Args:
        **kwargs: See `PyRosettaCluster` docstring. The keyword arguments must also include
            'protocols', an iterable object of function or generator objects specifying
            an ordered sequence of user-defined PyRosetta protocols to execute for
            the simulation (see `PyRosettaCluster().distribute` docstring). The keyword arguments
            may also optionally include 'clients_indices' or 'resources' (see
            `PyRosettaCluster().distribute` docstring).
    """
    protocols = kwargs.pop("protocols", None)
    clients_indices = kwargs.pop("clients_indices", None)
    resources = kwargs.pop("resources", None)
    PyRosettaCluster(**kwargs).distribute(
        protocols=protocols,
        clients_indices=clients_indices,
        resources=resources,
    )

run: Callable[..., Optional[NoReturn]] = produce

@wraps(produce, assigned=("__doc__",), updated=())
def iterate(**kwargs: Any) -> Union[NoReturn, Generator[Tuple[PackedPose, Dict[Any, Any]], None, None]]:
    protocols = kwargs.pop("protocols", None)
    clients_indices = kwargs.pop("clients_indices", None)
    resources = kwargs.pop("resources", None)
    for result in PyRosettaCluster(**kwargs).generate(
        protocols=protocols,
        clients_indices=clients_indices,
        resources=resources,
    ):
        yield result

produce.__doc__ += """
    Returns:
        None
    """
iterate.__doc__ = iterate.__doc__.replace(
    "PyRosettaCluster().distribute", "PyRosettaCluster().generate"
) + """
    Yields:
        (PackedPose, dict) tuples from the most recently run user-provided PyRosetta protocol if
        `PyRosettaCluster(save_all=True)` otherwise from the final user-defined PyRosetta protocol.
    """
