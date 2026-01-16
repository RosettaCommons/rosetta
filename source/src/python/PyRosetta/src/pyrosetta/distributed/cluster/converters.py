# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import git
    import toolz
    from dask.distributed import Client, LocalCluster
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.converters' requires the "
        + "third-party packages 'dask', 'gitpython' and 'toolz' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/dask/\n"
        + "https://gitpython.readthedocs.io/en/stable/intro.html\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import collections
import logging
import os
import pyrosetta
import sys
import types
import warnings

from functools import singledispatch
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import (
    Any,
    Callable,
    Dict,
    Iterable,
    List,
    NoReturn,
    Optional,
    Sized,
    Tuple,
    TypeVar,
    Union,
)

from pyrosetta.distributed.cluster.config import (
    get_environment_cmd,
    get_environment_manager,
)
from pyrosetta.distributed.cluster.converter_tasks import (
    get_yml,
    is_bytes,
    is_dict,
    is_empty,
    is_packed,
    parse_input_packed_pose as _parse_input_packed_pose,
    to_int,
    to_iterable,
    to_packed,
    to_str,
    maybe_issue_environment_warnings as _maybe_issue_environment_warnings,
)
from pyrosetta.distributed.cluster.serialization import Serialization


S = TypeVar("S", bound=Serialization)


def _parse_filter_results(obj: Any) -> Union[bool, NoReturn]:
    """Parse the input `filter_results` attribute of PyRosettaCluster."""
    _issue_future_warning = True

    @singledispatch
    def converter(obj: Any) -> NoReturn:
        raise ValueError("'filter_results' must be of type `bool` or `NoneType`!")

    @converter.register(type(None))
    def _parse_none(obj: None) -> bool:
        if _issue_future_warning:
            warnings.warn(
                (
                    "As of PyRosettaCluster version 2.1.0, the 'filter_results' instance attribute "
                    "is enabled by default, which automatically filters empty `PackedPose` objects between "
                    "user-provided PyRosetta protocols to help reduce compute overhead. Please explicitly set "
                    "either `filter_results=True` (the currently enabled, new setting) or `filter_results=False` "
                    "(to revert to legacy behavior before version 2.1.0) to silence this notice. This notice "
                    "will disappear in a future version of PyRosettaCluster."
                ),
                FutureWarning,
                stacklevel=5,
            )
        return True

    @converter.register(bool)
    def _parse_bool(obj: bool) -> bool:
        return obj

    return converter(obj)


def _parse_decoy_ids(objs: Any) -> List[int]:
    """
    Normalize user-provided PyRosetta 'decoy_ids' to a `list` object containing `int` objects.
    """

    return to_iterable(objs, to_int, "decoy_ids")


def _parse_empty_queue(protocol_name: str, ignore_errors: bool) -> None:
    """Return a `None` object when a protocol results in an error with `ignore_errors=True`."""
    logging.warning(
        f"User-provided PyRosetta protocol '{protocol_name}' resulted in an empty queue with `ignore_errors={ignore_errors}`! "
        + "Putting a `None` object into the queue."
    )
    return None


def _parse_environment(obj: Any) -> Union[str, NoReturn]:
    """Parse the input `environment` attribute of PyRosettaCluster."""

    @singledispatch
    def converter(obj: Any) -> NoReturn:
        raise ValueError("The 'environment' instance attribute must be of type `str` or `NoneType`!")

    @converter.register(type(None))
    def _parse_none(obj: None) -> str:
        yml = get_yml()
        if yml == "":
            environment_cmd = get_environment_cmd()
            environment_manager = get_environment_manager()
            _warning_msg = (
                "`{0}` did not run successfully, "
                "so the active {1} file string was not saved! "
                "It is recommended to run:\n`{2}`\n"
                "to reproduce this simulation later."
            )
            if environment_manager == "pixi":
                logging.warning(
                    _warning_msg.format(
                        environment_cmd,
                        "pixi project lock",
                        environment_cmd,
                    )
                )
            elif environment_manager == "uv":
                logging.warning(
                    _warning_msg.format(
                        environment_cmd,
                        "uv project requirements",
                        f"{environment_cmd} > requirements.txt",
                    )
                )
            elif environment_manager == "mamba":
                logging.warning(
                    _warning_msg.format(
                        environment_cmd,
                        "mamba environment YML",
                        f"{environment_cmd} > environment.yml",
                    )
                )
            elif environment_manager == "conda":
                logging.warning(
                    _warning_msg.format(
                        environment_cmd,
                        "conda environment YML",
                        f"{environment_cmd} > environment.yml",
                    )
                )
            else:
                raise RuntimeError(f"Unsupported environment manager: {environment_manager}")
        return yml

    @converter.register(str)
    def _parse_str(obj: str) -> Union[str, NoReturn]:
        environment_manager = get_environment_manager()
        if obj == "":
            _warning_msg = (
                "The input 'environment' parameter argument is an empty string, "
                "which is not a valid {0} file string capturing the active {1}! "
                "Reproduction simulations may not necessarily reproduce "
                "the original decoy(s)! Please verify that your active "
                "{1} is identical to the original {1} that "
                "generated the decoy(s) you wish to reproduce!"
                "\nBypassing {1} validation...\n"
            )
            if environment_manager == "pixi":
                logging.warning(_warning_msg.format("pixi lock", "pixi project"))
            elif environment_manager == "uv":
                logging.warning(_warning_msg.format("uv requirements", "uv project"))
            elif environment_manager == "mamba":
                logging.warning(_warning_msg.format("YML", "mamba environment"))
            elif environment_manager == "conda":
                logging.warning(_warning_msg.format("YML", "conda environment"))
            else:
                raise RuntimeError(f"Unsupported environment manager: {environment_manager}")
            return obj
        else:
            if obj != get_yml():
                _err_msg = (
                    "The 'environment' parameter argument is not equivalent to the {0} file string "
                    "generated by the active {1}, and therefore the original "
                    "decoy may not necessarily be reproduced. Please set the 'environment' parameter "
                    "argument to an empty string ('') to bypass {1} validation and run the simulation."
                )
                if environment_manager == "pixi":
                    raise AssertionError(_err_msg.format("pixi lock", "pixi project"))
                elif environment_manager == "uv":
                    raise AssertionError(_err_msg.format("uv requirements", "uv project"))
                elif environment_manager == "mamba":
                    raise AssertionError(_err_msg.format("YML", "mamba environment"))
                elif environment_manager == "conda":
                    raise AssertionError(_err_msg.format("YML", "conda environment"))
                else:
                    raise RuntimeError(f"Unsupported environment manager: {environment_manager}")
            else:
                _debug_msg = "The 'environment' parameter argument correctly validated against the active {0}!"
                if environment_manager in ("pixi", "uv"):
                    logging.debug(_debug_msg.format(f"{environment_manager} project"))
                else:
                    logging.debug(_debug_msg.format(f"{environment_manager} environment"))
                return obj

    return converter(obj)


def _parse_protocols(objs: Any) -> Union[List[Union[Callable[..., Any], Iterable[Any]]], NoReturn]:
    """
    Parse the `protocols` argument parameters from the PyRosettaCluster().distribute() method.
    """

    @singledispatch
    def converter(objs: Any) -> NoReturn:
        raise RuntimeError(
            "The user-provided PyRosetta protocols must be an "
            + "iterable of objects of `types.GeneratorType` and/or "
            + "`types.FunctionType` types, or an instance of type "
            + "`types.GeneratorType` or `types.FunctionType`, not of type "
            + f"{type(objs)}!"
        )

    @converter.register(type(None))
    def _none_to_list(obj: None) -> List[Any]:
        return []

    @converter.register(types.FunctionType)
    @converter.register(types.GeneratorType)
    def _func_to_list(
        obj: Union[Callable[..., Any], Iterable[Any]]
    ) -> List[Union[Callable[..., Any], Iterable[Any]]]:
        return [obj]

    @converter.register(collections.abc.Iterable)
    def _to_list(
        objs: Iterable[Any],
    ) -> Union[List[Union[Callable[..., Any], Iterable[Any]]], NoReturn]:
        for obj in objs:
            if not isinstance(obj, (types.FunctionType, types.GeneratorType)):
                raise TypeError(
                    "Each member of PyRosetta protocols must be of type "
                    + "`types.FunctionType` or `types.GeneratorType`! "
                    + f"Received: {type(obj)}"
                )
        return list(objs)

    return converter(objs)


def _parse_yield_results(yield_results: Any) -> Union[bool, NoReturn]:
    @singledispatch
    def converter(objs: Any) -> NoReturn:
        raise ValueError("'yield_results' parameter must be of type `bool`.")
    
    @converter.register(bool)
    def _to_bool(objs: bool) -> bool:
        return objs
    
    @converter.register(int)
    @converter.register(float)
    @converter.register(str)
    def _int_to_bool(objs: Union[int, float, str]) -> bool:
        return bool(objs)
    
    @converter.register(type(None))
    def _none_to_bool(objs: None) -> bool:
        return False
    
    return converter(yield_results)


def _parse_norm_task_options(obj: Any) -> Union[bool, NoReturn]:
    """Parse the input `norm_task_options` attribute of PyRosettaCluster."""
    _issue_future_warning = True

    @singledispatch
    def converter(obj: Any) -> NoReturn:
        raise ValueError("'norm_task_options' must be of type `bool` or `NoneType`!")

    @converter.register(type(None))
    def _parse_none(obj: None) -> bool:
        if _issue_future_warning:
            warnings.warn(
                (
                    "As of PyRosettaCluster version 2.3.0, the 'norm_task_options' instance attribute is enabled by "
                    "default, which automatically normalizes the task 'options' and 'extra_options' keyword arguments "
                    "for facile reproducibility. Please explicitly set either `norm_task_options=True` (the currently "
                    "enabled, new setting) or `norm_task_options=False` (to revert to legacy behavior before version 2.3.0) "
                    "to silence this notice. This notice will disappear in a future version of PyRosettaCluster."
                ),
                FutureWarning,
                stacklevel=5,
            )
        return True

    @converter.register(bool)
    def _parse_bool(obj: bool) -> bool:
        return obj

    return converter(obj)


def _parse_pyrosetta_build(obj: Any) -> Union[str, NoReturn]:
    """Parse the PyRosetta build string."""

    _pyrosetta_version_string = pyrosetta._build_signature()

    @singledispatch
    def converter(obj: Any) -> NoReturn:
        raise ValueError("'pyrosetta_build' must be of type `str` or `NoneType`!")

    @converter.register(type(None))
    def _default_none(obj: None) -> str:
        return _pyrosetta_version_string

    @converter.register(str)
    def _validate_pyrosetta_version_string(obj: str) -> Union[str, NoReturn]:
        if obj == "":
            logging.warning(
                "The input 'pyrosetta_build' keyword argument parameter is an empty string, "
                + "which is not a valid PyRosetta build string capturing the original PyRosetta "
                + "version! Reproduction simulations may not necessarily reproduce "
                + "the original decoy(s)! Please verify that your PyRosetta build "
                + "is identical to the original PyRosetta build that "
                + "generated the decoy(s) you wish to reproduce!"
                + "\nBypassing PyRosetta build validation...\n"
            )
            return obj
        else:
            if obj != _pyrosetta_version_string:
                raise AssertionError(
                    f"The user-provided PyRosetta build string '{obj}' is not equal "
                    + f"to the currently used PyRosetta build string '{_pyrosetta_version_string}'! "
                    + "Therefore, the original decoy may not necessarily be reproduced! "
                    + "Using different PyRosetta builds can lead to different outputs. "
                    + "Please consider running this simulation using the PyRosetta build that "
                    + "was used with the original simulation run. Please set the 'pyrosetta_build' keyword "
                    + "argument parameter to an empty string ('') to bypass PyRosetta build validation."
                )
            else:
                logging.debug(
                    "The 'pyrosetta_build' keyword argument parameter correctly validated against the original PyRosetta build!"
                )
                return obj

    return converter(obj)


def _parse_scratch_dir(obj: Any) -> Union[str, NoReturn]:
    """Parse the input `scratch_dir` attribute of PyRosettaCluster."""

    @singledispatch
    def converter(obj: Any) -> NoReturn:
        raise ValueError(
            f"'scratch_dir' directory {obj} could not be created or found!"
        )

    @converter.register(str)
    def _is_str(obj: str) -> str:
        return obj

    @converter.register(type(None))
    def _default_none(obj: None) -> str:
        temp_path = os.sep + "temp"
        if os.path.exists(temp_path):
            scratch_dir = temp_path
        else:
            scratch_dir = os.getcwd()
        return scratch_dir

    return converter(obj)


def _parse_seeds(objs: Any) -> List[str]:
    """
    Normalize user-provided PyRosetta 'seeds' to a `list` object containing `str` objects.
    """

    return to_iterable(objs, to_str, "seeds")


def _parse_sha1(obj: Any) -> Union[str, NoReturn]:
    """Parse the input `sha1` attribute of PyRosettaCluster."""

    @singledispatch
    def converter(obj: Any) -> NoReturn:
        raise ValueError(f"'sha1' attribute must be of type `str` or `NoneType`!")

    @converter.register(str)
    def _register_str(obj: str) -> Union[str, NoReturn]:
        """Parse `str` type inputs to the 'sha1' attribute."""

        try:
            repo = git.Repo(".", search_parent_directories=True)
        except git.InvalidGitRepositoryError as ex:
            raise git.InvalidGitRepositoryError(
                f"{ex}\nThe script being executed is not in a git repository! "
                + "It is strongly recommended to use version control for your "
                + "scripts. To continue without using version control, set the "
                + "`sha1` attribute of PyRosettaCluster to `None`."
            )
        # Path to and name of the git repository
        git_root = repo.git.rev_parse("--show-toplevel")
        repo_name = os.path.split(git_root)[-1]
        try:
            commit = repo.head.commit
        except ValueError as ex:
            raise git.InvalidGitRepositoryError(
                f"{ex}\nNo HEAD commit present! Is this a repository with no commits?"
            )
        if obj == "":
            # Ensure that everything in the work directory is included in the repository
            if len(repo.untracked_files):
                raise git.InvalidGitRepositoryError(
                    "There are untracked files in your git repository! "
                    + "If these are important for the simulation, you will not "
                    + "be able to reproduce your work! Commit the changes if "
                    + "appropriate, or add untracked files to your .gitignore file "
                    + "before running PyRosettaCluster!"
                )
            # Ensure that all changes to tracked files are committed
            if repo.is_dirty(untracked_files=True):
                raise git.InvalidGitRepositoryError(
                    "The working directory is dirty! "
                    + "Commit local changes to ensure reproducibility."
                )
            return commit.hexsha
        else:
            # A sha1 was provided. Validate that it is HEAD
            if not commit.hexsha.startswith(obj):
                logging.error(
                    "A `sha1` attribute was provided to PyRosettaCluster, "
                    + "but it appears that you have not checked out this revision!"
                    + f"Please run on the command line:\n`git checkout {obj}`\n"
                    + "and then re-execute this script."
                )
                raise RuntimeError(
                    "The `sha1` attribute provided to PyRosettaCluster "
                    + "does not match the current `HEAD`! "
                    + "See log files for details on resolving the issue."
                )
            return obj

    @converter.register(type(None))
    def _default_none(obj: None) -> str:
        """If the user provides `None` to the `sha1` attribute, return an empty string."""

        return ""

    return converter(obj)


def _parse_system_info(obj: Any) -> Union[Dict[Any, Any], NoReturn]:
    """Parse the input `system_info` attribute of PyRosettaCluster."""

    _sys_platform = {"sys.platform": sys.platform}

    @singledispatch
    def converter(obj: Any) -> NoReturn:
        raise ValueError("'system_info' must be of type `dict` or `NoneType`!")

    @converter.register(type(None))
    def _default_none(obj: None) -> Dict[str, str]:
        return _sys_platform

    @converter.register(dict)
    def _overwrite_sys_platform(obj: Dict[Any, Any]) -> Dict[Any, Any]:
        if "sys.platform" in obj:
            _obj_sys_platform = obj["sys.platform"]
            _curr_sys_platform = _sys_platform["sys.platform"]
            if _obj_sys_platform != _curr_sys_platform:
                logging.warning(
                    "The dictionary key 'sys.platform' of the PyRosettaCluster "
                    + "'system_info' attribute indicates a previously used system "
                    + f"platform '{_obj_sys_platform}', but the currently used "
                    + f"system platform is '{_curr_sys_platform}'! Therefore, the original "
                    + "decoy may not necessarily be reproduced! Platform-specific "
                    + "information can lead to different outputs. Please consider "
                    + "running this simulation in a Docker container with the same "
                    + "system platform as the original simulation run, or switching "
                    + "to the original system platform before reproducing this simulation."
                )
        return toolz.dicttoolz.merge(obj, _sys_platform)

    return converter(obj)


def _parse_logging_address(self) -> str:
    _default_local = "localhost:0"
    _default_remote = "0.0.0.0:0"

    if all(attribute is None for attribute in (self.scheduler, self.client, self.clients)):
        logging_address = _default_local
    elif isinstance(self.client, Client):
        logging_address = (
            _default_local if isinstance(self.client.cluster, LocalCluster) else _default_remote
        )
    elif all(isinstance(_client, Client) for _client in self.clients):
        logging_address = (
            _default_local
            if all(isinstance(_client.cluster, LocalCluster) for _client in self.clients)
            else _default_remote
        )
    else:
        logging_address = _default_remote

    return logging_address


def _get_decoy_id(protocols: Sized, decoy_ids: List[int]) -> Optional[int]:
    """Get the decoy number given the user-provided PyRosetta protocols."""

    if decoy_ids:
        decoy_id_index = (len(decoy_ids) - len(protocols)) - 1
        decoy_id = decoy_ids[decoy_id_index]
    else:
        decoy_id = None

    return decoy_id


def _get_packed_poses_output_kwargs(
    result: Any,
    input_kwargs: Dict[Any, Any],
    protocol_name: str,
) -> Tuple[List[PackedPose], Dict[Any, Any]]:
    packed_poses = []
    protocol_kwargs = []
    for obj in to_iterable(result, to_packed, protocol_name):
        if is_packed(obj):
            packed_poses.append(obj)
        elif is_dict(obj):
            protocol_kwargs.append(obj)

    if len(packed_poses) == 0:
        packed_poses = to_iterable(None, to_packed, protocol_name)

    if len(protocol_kwargs) == 0:
        output_kwargs = input_kwargs
    elif len(protocol_kwargs) == 1:
        output_kwargs = next(iter(protocol_kwargs))
        output_kwargs.update(
            toolz.dicttoolz.keyfilter(lambda k: k.startswith("PyRosettaCluster_"), input_kwargs)
        )
    elif len(protocol_kwargs) >= 2:
        raise ValueError(
            f"User-provided PyRosetta protocol '{protocol_name}' may return at most one object of type `dict`."
        )

    return packed_poses, output_kwargs


def _get_compressed_packed_pose_kwargs_pairs_list(
    packed_poses: List[PackedPose],
    output_kwargs: Dict[Any, Any],
    protocol_name: str,
    protocols_key: str,
    decoy_ids: List[int],
    serializer: S,
) -> List[Tuple[bytes, bytes]]:
    decoy_id = _get_decoy_id(output_kwargs[protocols_key], decoy_ids)
    compressed_packed_pose_kwargs_pairs_list = []
    for i, packed_pose in enumerate(packed_poses):
        if (decoy_id != None) and (i != decoy_id):
            logging.info(
                "Discarding a returned decoy because it does not match the user-provided 'decoy_ids'."
            )
            continue
        task_kwargs = serializer.deepcopy_kwargs(output_kwargs)
        if "PyRosettaCluster_decoy_ids" not in task_kwargs:
            task_kwargs["PyRosettaCluster_decoy_ids"] = []
        task_kwargs["PyRosettaCluster_decoy_ids"].append((protocol_name, i))
        compressed_packed_pose = serializer.compress_packed_pose(packed_pose)
        compressed_task_kwargs = serializer.compress_kwargs(task_kwargs)
        compressed_packed_pose_kwargs_pairs_list.append((compressed_packed_pose, compressed_task_kwargs))

    if decoy_id:
        assert (
            len(compressed_packed_pose_kwargs_pairs_list) == 1
        ), "When specifying decoy_ids, there may only be one decoy_id per protocol."

    return compressed_packed_pose_kwargs_pairs_list


def _parse_protocol_results(
    result: Any,
    input_kwargs: Dict[Any, Any],
    protocol_name: str,
    protocols_key: str,
    decoy_ids: List[int],
    serializer: S,
) -> List[Tuple[bytes, bytes]]:
    """Parse results from the user-provided PyRosetta protocol."""
    packed_poses, output_kwargs = _get_packed_poses_output_kwargs(result, input_kwargs, protocol_name)
    compressed_packed_pose_kwargs_pairs_list = _get_compressed_packed_pose_kwargs_pairs_list(
        packed_poses, output_kwargs, protocol_name, protocols_key, decoy_ids, serializer
    )

    return compressed_packed_pose_kwargs_pairs_list


def _parse_target_results(objs: List[Tuple[bytes, bytes]]) -> List[Tuple[bytes, bytes]]:
    """Parse results returned from the spawned thread."""

    ids = set()
    n_obj = 0
    for obj in objs:
        if len(obj) != 2:
            raise TypeError("Returned result should be of length 2.")
        n_obj += len(obj)
        ids.update(set(map(id, obj)))
    assert len(ids) == n_obj, "Returned results do not have unique memory addresses."

    return objs


def _parse_tasks(objs: Any) -> Union[List[Dict[Any, Any]], NoReturn]:
    """Parse the input `tasks` attribute of PyRosettaCluster."""

    @singledispatch
    def converter(objs: Any) -> NoReturn:
        raise ValueError(
            "Parameter passed to 'tasks' argument must be an iterable, a function, a generator, or None!"
        )

    @converter.register(dict)
    def _from_dict(obj: Dict[Any, Any]) -> List[Dict[Any, Any]]:
        return [obj]

    @converter.register(type(None))
    def _from_none(obj: None) -> List[Dict[Any, Any]]:
        logging.warning(
            "PyRosettaCluster `tasks` attribute was set to `None`! Using a default empty task."
        )
        return [{}]

    @converter.register(types.FunctionType)
    def _from_function(
        obj: Callable[..., Iterable[Any]]
    ) -> Union[List[Dict[Any, Any]], NoReturn]:
        _tasks = []
        for obj in objs():
            if isinstance(obj, dict):
                _tasks.append(obj)
            else:
                raise TypeError(
                    f"Each task must be of type `dict`. Received '{obj}' of type `{type(obj)}`."
                )
        return _tasks

    @converter.register(collections.abc.Iterable)
    @converter.register(types.GeneratorType)
    def _from_iterable(obj: Iterable[Any]) -> Union[List[Dict[Any, Any]], NoReturn]:
        _tasks = []
        for obj in objs:
            if isinstance(obj, dict):
                _tasks.append(obj)
            else:
                raise TypeError(
                    f"Each task must be of type `dict`. Received '{obj}' of type `{type(obj)}`."
                )
        return _tasks

    return converter(objs)


def _parse_output_decoy_types(objs: Any) -> Union[List[str], NoReturn]:
    """Parse the input `output_decoy_types` attribute of PyRosettaCluster."""
    _output_decoy_types = (".pdb", ".pkl_pose", ".b64_pose", ".init")

    @singledispatch
    def converter(objs: Any) -> NoReturn:
        raise ValueError(
            "Parameter passed to 'output_decoy_types' argument must be "
            + "an iterable of `str` objects, or `None`."
        )

    @converter.register(type(None))
    def _from_none(obj: None) -> List[str]:
        return [_output_decoy_types[0]]

    @converter.register(collections.abc.Iterable)
    def _from_iterable(objs: collections.abc.Iterable) -> Union[List[str], NoReturn]:
        _seen = set()
        _types = []
        for obj in objs:
            if isinstance(obj, str):
                if obj in _output_decoy_types and obj not in _seen:
                    _types.append(obj)
                    _seen.add(obj)
                else:
                    raise ValueError(f"Available output decoy types are: {_output_decoy_types}")
            else:
                converter.dispatch(object)(obj)

        return _types if len(_types) >= 1 else converter(None)

    return converter(objs)


def _parse_output_scorefile_types(objs: Any) -> Union[List[str], NoReturn]:
    """Parse the input `output_scorefile_types` attribute of PyRosettaCluster."""
    _output_scorefile_types = (".json",)

    @singledispatch
    def converter(objs: Any) -> NoReturn:
        raise ValueError(
            "Parameter passed to 'output_scorefile_types' argument must be "
            + "an iterable of `str` objects, or `None`."
        )

    @converter.register(type(None))
    def _from_none(obj: None) -> List[str]:
        return [_output_scorefile_types[0]]

    @converter.register(collections.abc.Iterable)
    def _from_iterable(objs: collections.abc.Iterable) -> Union[List[str], NoReturn]:
        _seen = set()
        _types = []
        for obj in objs:
            if isinstance(obj, str) and obj not in _seen:
                if not obj.startswith("."):
                    raise ValueError(
                        f"The element '{obj}' in the 'output_scorefile_types' keyword argument "
                        "parameter must start with '.' to be a valid filename extension."
                    )
                _types.append(obj)
                _seen.add(obj)
            else:
                converter.dispatch(object)(obj)

        return _types if len(_types) >= 1 else converter(None)

    result = converter(objs)
    if result != converter(None) and "pandas" not in pyrosetta.secure_unpickle.get_secure_packages():
        _compressions = [f"'{t}'" for t in result if t != _output_scorefile_types[0]]
        if len(_compressions) > 1:
            _compressions[-1] = f"and {_compressions[-1]}"
        _compression_msg = ", ".join(_compressions) if len(_compressions) > 2 else " ".join(_compressions)
        _msg = (
            f"In order to save PyRosettaCluster scorefiles with {_compression_msg} file type compression, "
            + "the output `pandas.DataFrame()` objects must be securely unpicklable. In order to securely "
            + "unpickle `pandas.DataFrame()` objects in PyRosetta, please add 'pandas' as a secure package "
            + "using `pyrosetta.secure_unpickle.add_secure_package('pandas')`, then try again!"
        )
        raise AssertionError(_msg)

    return result


def _version_tuple_to_str(version_tuple: Tuple[int, int, int]) -> str:
    """Return a version string from a version tuple."""
    return ".".join(map(str, version_tuple))
