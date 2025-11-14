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
        "Importing 'pyrosetta.distributed.cluster.converter_tasks' requires the "
        + "third-party packages 'distributed' and 'toolz' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/distributed/\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import bz2
import collections
import json
import logging
import os
import pyrosetta.distributed.io as io
import re
import shutil
import subprocess
import warnings

from contextlib import contextmanager
from functools import singledispatch
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.exceptions import PyRosettaIsNotInitializedError
from pyrosetta.rosetta.basic import was_init_called
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.utility.initialization import PyRosettaInitDictWriter
from typing import (
    Any,
    Callable,
    Dict,
    Generator,
    Iterable,
    List,
    NoReturn,
    Optional,
    Tuple,
    TypeVar,
    Union,
)
from urllib.parse import urlparse, urlunparse

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
from pyrosetta.distributed.cluster.io import (
    IO,
    get_poses_from_init_file,
    secure_read_pickle,
    sign_init_file_metadata_and_poses,
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
            yml = get_yml()
            if yml == "":
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
            elif "pyrosetta=" not in yml:
                warnings.warn(
                    "The currently installed 'pyrosetta' package version is not specified in the exported environment file! "
                    + "Consequently, the PyRosettaCluster simulation will be difficult to reproduce at a later time. "
                    + "To use the `pyrosetta.distributed.cluster` namespace and ensure reproducibility of PyRosetta simulations, "
                    + "please re-install the 'pyrosetta' package using the Rosetta Commons conda channel. For instructions, visit:\n"
                    + "https://www.pyrosetta.org/downloads\nNote that installing PyRosetta using pip and the 'pyrosetta-installer' "
                    + "package does not pin the PyRosetta version to the currently activated virtual environment.",
                    UserWarning,
                    stacklevel=4,
                )  # Warn that the PyRosetta package version is not specified in the active virtual environment
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
    input_file: Optional[Union[str, Pose, PackedPose]] = None,
    scorefile: Optional[str] = None,
    decoy_name: Optional[str] = None,
) -> Union[List[str], NoReturn]:
    """
    Get the user-defined PyRosetta protocols as a `list` object of `str` objects.

    Args:
        input_file: A `str` object specifying the path to the '.pdb', '.pdb.bz2', '.pkl_pose',
            '.pkl_pose.bz2', '.b64_pose', '.b64_pose.bz2', '.init', or '.init.bz2' file,
            or a `Pose` or `PackedPose` object, from which to extract PyRosettaCluster instance
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
            warnings.warn(
                "Received 'input_file' and either 'scorefile' or 'decoy_name' keyword argument parameters. "
                + "Ignoring 'scorefile' and 'decoy_name' and using 'input_file' keyword argument parameter!",
                UserWarning,
                stacklevel=3,
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


def get_scores_dict(obj: Union[str, Pose, PackedPose]) -> Union[Dict[str, Dict[str, Any]], NoReturn]:
    """
    Get the PyRosettaCluster scores dictionary from either a `Pose` or `PackedPose` object, or a '.pdb',
    '.pdb.bz2', '.pkl_pose', '.pkl_pose.bz2', '.b64_pose', '.b64_pose.bz2', '.init', or '.init.bz2' file.
    """

    if isinstance(obj, (Pose, PackedPose)):
        pdbstring = io.to_pdbstring(obj)
    elif isinstance(obj, str):
        if not os.path.exists(obj):
            raise IOError(f"The `input_file` argument parameter must exist on disk! Received: '{obj}'")
        if obj.endswith(".pdb.bz2"):
            with open(obj, "rb") as fbz2:
                pdbstring = bz2.decompress(fbz2.read()).decode()
        elif obj.endswith(".pdb"):
            with open(obj, "r") as f:
                pdbstring = f.read()
        elif obj.endswith((".pkl_pose", ".pkl_pose.bz2", ".b64_pose", ".b64_pose.bz2")):
            if not was_init_called():
                raise PyRosettaIsNotInitializedError(
                    "To get the PyRosettaCluster scores dictionary from a '.pkl_pose', '.pkl_pose.bz2', "
                    + "'.b64_pose' or '.b64_pose.bz2' file, PyRosetta must be initialized (with the same "
                    + "residue type set that was used to save the original decoy output file)."
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
        elif obj.endswith((".init", ".init.bz2")):
            if not was_init_called():
                if obj.endswith(".init.bz2"):
                    raise PyRosettaIsNotInitializedError(
                        "To get the PyRosettaCluster scores dictionary from a '.init.bz2' file, please first initialize "
                        + f"PyRosetta using the `pyrosetta.distributed.io.init_from_file()` function: '{obj}'"
                    )
                elif obj.endswith(".init"):
                    raise PyRosettaIsNotInitializedError(
                        "To get the PyRosettaCluster scores dictionary from a '.init' file, please first initialize "
                        + f"PyRosetta using the `pyrosetta.init_from_file()` function: '{obj}'"
                    )
            _input_packed_pose, output_packed_pose = get_poses_from_init_file(obj, verify=True)
            if output_packed_pose is None:
                raise ValueError(
                    "The input '.init' or '.init.bz2' file does not contain an output decoy from a PyRosettaCluster simulation: "
                    + f"'{obj}'. To get the PyRosettaCluster scores dictionary from a '.init' or '.init.bz2' file, please ensure "
                    + "that `pyrosetta.distributed.cluster.export_init_file()` was run on the original decoy output file, "
                    + "or that the `PyRosettaCluster(output_decoy_types=['.init'])` PyRosetta initialization file output "
                    + "decoy type was enabled in the original PyRosettaCluster simulation."
                )
            pdbstring = io.to_pdbstring(output_packed_pose)
        else:
            raise ValueError(
                "The `input_file` argument parameter must end in '.pdb', '.pdb.bz2', '.pkl_pose', '.pkl_pose.bz2', "
                + "'.b64_pose', '.b64_pose.bz2', '.init', or '.init.bz2'. "
                + f"Received: '{obj}'"
            )
    else:
        raise TypeError(
            "The `input_file` argument parameter must be a `Pose`, `PackedPose` or `str` object that ends in '.pdb', "
            + "'.pdb.bz2', '.pkl_pose', '.pkl_pose.bz2', '.b64_pose', '.b64_pose.bz2', '.init', or '.init.bz2'. "
            + f"Received: '{type(obj)}'"
        )

    scores_dict = None
    for line in reversed(pdbstring.split(os.linesep)):
        if line.startswith(IO.REMARK_FORMAT):
            scores_dict = json.loads(
                line.split(IO.REMARK_FORMAT)[-1],
                cls=None,
                object_hook=None,
                object_pairs_hook=None,
                parse_int=None,
                parse_constant=None,
            )
            break
    else:
        _err_msg = f"Could not parse '{IO.REMARK_FORMAT}' comment from the input object: '{obj}'"
        if isinstance(obj, (Pose, PackedPose)):
            raise ValueError(
                f"{_err_msg}. If the '{type(obj)}' object was initialized from a '.pdb' or '.pdb.bz2' "
                + "file output by PyRosettaCluster, please input the file path directly into the "
                + f"`get_scores_dict()` function to parse the '{IO.REMARK_FORMAT}' comment."
            )
        else:
            raise IOError(_err_msg)

    if scores_dict is None:
        raise IOError(f"Could not parse the input argument parameter: '{obj}'")
    if not all(d in scores_dict for d in ("instance", "metadata", "scores")):
        raise KeyError(f"Could not parse the input argument parameter: '{obj}'")

    return scores_dict


def export_init_file(
    output_file: str,
    output_init_file: Optional[str] = None,
    compressed: Optional[bool] = None,
) -> None:
    """
    Export a PyRosetta initialization file from a decoy output file. The PyRosettaCluster simulation
    that produced the decoy output file must have had the 'output_init_file' instance attribute enabled,
    so the 'init_file' key value can be detected in the metadata of the decoy output file. This function
    can be used to prepend the decoy output file to the detected PyRosetta initialization file for more
    facile simulation reproduction using the `reproduce()` function.

    Args:
        output_file: A required `str` object representing the decoy output file. The file must end in
            either: '.pdb', '.pdb.bz2', '.pkl_pose', '.pkl_pose.bz2', '.b64_pose', or '.b64_pose.bz2'.
        output_init_file: An optional `str` object specifying the output PyRosetta initialization file
            path ending with '.init'. If `NoneType` is provided, then the PyRosetta initialization file
            path is derived from the 'output_file' argument parameter by replacing the file extension
            with '.init' (or '.init.bz2' when the 'compressed' argument parameter is set to `True`).
            Default: None
        compressed: A `bool` object specifying whether or not to compress the output PyRosetta initialization
            file with `bzip2`, resulting in a '.init.bz2' output PyRosetta initialization file.
                Default: True

    Returns:
        None
    """

    _types = (".pdb", ".pdb.bz2", ".pkl_pose", ".pkl_pose.bz2", ".b64_pose", ".b64_pose.bz2")

    if isinstance(output_file, str) and os.path.isfile(output_file) and output_file.endswith(_types):
        if output_init_file is None:
            for _type in _types:
                if output_file.endswith(_type):
                    output_init_file = f"{output_file[: -len(_type)]}.init"
                    break
        elif isinstance(output_init_file, str):
            if not output_init_file.endswith(".init"):
                raise ValueError(
                    "The 'output_init_file' keyword argument parameter must end with '.init'. "
                    f"Received: '{output_init_file}'"
            )
        else:
            raise TypeError(
                "The 'output_init_file' keyword argument parameter must be a `str` or `NoneType` object. "
                + f"Received: {type(output_init_file)}"
            )
        if not isinstance(compressed, (type(None), bool)):
            raise TypeError(
                "The 'compressed' keyword argument parameter must be a `bool` or `NoneType` object. "
                + f"Received: {type(compressed)}"
            )
        if compressed:
            output_init_file += ".bz2"
        if os.path.isfile(output_init_file):
            raise FileExistsError(
                f"The PyRosetta initialization file path already exists: '{output_init_file}'. "
                + "Please set the 'output_init_file' keyword argument parameter to a different value."
            )

        scores_dict = get_scores_dict(output_file)
        init_file = scores_dict["metadata"]["init_file"]
        if init_file:
            if os.path.isfile(init_file):
                if not was_init_called():
                    raise PyRosettaIsNotInitializedError(
                        "In order to export a PyRosetta initialization file, please ensure that PyRosetta is already "
                        + "initialized (using the `pyrosetta.distributed.io.init_from_file()` function) with the "
                        + f"following PyRosetta initialization file, and then run `export_init_file()`: '{init_file}'"
                    )
                input_packed_pose, _output_packed_pose = get_poses_from_init_file(init_file, verify=True)
                if _output_packed_pose is not None:
                    raise ValueError(
                        "The 'output_file' argument parameter already contains an output decoy in the "
                        + f"detected PyRosetta initialization file: '{init_file}'. Aborting export!"
                    )
                if output_file.endswith(".bz2"):
                    with open(output_file, "rb") as fbz2:
                        string = bz2.decompress(fbz2.read()).decode()
                    if output_file.endswith(".pdb.bz2"):
                        output_packed_pose = io.pose_from_pdbstring(string)
                    else:
                        output_packed_pose = io.to_packed(io.to_pose(string))
                else:
                    output_packed_pose = io.pose_from_file(output_file)
                # Cache simulation data from '.pdb' and '.pdb.bz2' files
                # loaded from `io.pose_from_pdbstring` or `io.pose_from_file`
                if output_file.endswith((".pdb", ".pdb.bz2")):
                    output_packed_pose = IO._add_pose_comment(
                        output_packed_pose,
                        IO._dump_json(get_scores_dict(output_file)),
                    )
                # Setup metadata and poses
                metadata, poses = sign_init_file_metadata_and_poses(
                    input_packed_pose=input_packed_pose,
                    output_packed_pose=output_packed_pose,
                )
                # Update PyRosetta initialization file
                init_dict = io.read_init_file(init_file)
                init_dict["metadata"] = metadata  # Overwrite
                init_dict["poses"] = poses  # Overwrite
                init_dict.update(dict(dry_run=False, overwrite=False, verbose=True))
                writer = PyRosettaInitDictWriter(**init_dict)
                init_file_json = writer.get_json() # Sign MD5
                if init_dict["verbose"]:
                    writer.print_cached_files(output_init_file, init_dict["dry_run"])
                if compressed:
                    with open(output_init_file, "wb") as f:
                        f.write(bz2.compress(str.encode(init_file_json)))
                else:
                    with open(output_init_file, "w") as f:
                        f.write(init_file_json)
                print(
                    f"Exported PyRosettaCluster decoy output file '{output_file}' to "
                    + f"PyRosetta initialization file: '{output_init_file}'"
                )
            else:
                raise ValueError(
                    "The 'output_file' argument parameter contains an 'init_file' key value in the "
                    + "cached metadata, but the specified '.init' or '.init.bz2' file cannot be found: "
                    + f"'{init_file}'. Please ensure the '.init' or '.init.bz2' file exists in the same "
                    + "path to use the `export_init_file()` function."
                )
        else:
            raise ValueError(
                "The 'output_file' argument parameter does not contain an 'init_file' key value "
                + "in the cached metadata, so the original simulation disabled the output of a "
                + "'.init' file. Exporting a '.init' file containing the 'output_file' argument "
                + "parameter is not supported without saving the original '.init' file. Please "
                + "enable the 'output_init_file' instance attribute in future PyRosettaCluster "
                + "simulations to use the `export_init_file()` function."
            )
    else:
        raise ValueError(
            "The 'output_file' argument parameter must be a `str` object, must exist on disk, and must "
            + f"end with one of the following filetype extensions: {_types}. Received: '{output_file}'"
        )


def get_yml() -> str:
    """
    Export the current environment to a string depending on the environment manager.
    """

    def remove_comments(text: str) -> str:
        """Remove lines starting with '#'."""
        return "\n".join(
            line for line in text.splitlines() if not line.strip().startswith("#")
        )

    def remove_metadata(text: str) -> str:
        """Remove 'name:' and 'prefix:' lines."""
        filtered_lines = [
            line
            for line in text.splitlines()
            if not line.startswith(("name:", "prefix:")) and line.strip()
        ]
        return "\n".join(filtered_lines) + "\n"

    def sanitize_url(url: str) -> str:
        """Remove username and password from URLs pointing to source domains."""
        parsed = urlparse(url)

        # No credentials present
        if "@" not in parsed.netloc:
            return url

        # Split credentials from host
        _credentials, host = parsed.netloc.split("@", 1)
        host_domain = host.split(":", 1)[0]  # Remove port if present

        # Only sanitize if domain matches sensitive domains
        if host_domain not in source_domains:
            return url

        # Build sanitized URL
        sanitized = parsed._replace(netloc=host)
        sanitized_url = urlunparse(sanitized)

        # Warn without leaking credentials
        warnings.warn(
            (
                "PyRosettaCluster automatically removed embedded credentials from the "
                f"conda channel '{host_domain}' while processing the environment file. "
                "These credentials are no longer required by this conda channel. "
                "Please remove them from your configuration to silence this warning."
            ),
            UserWarning,
            stacklevel=2,
        )

        return sanitized_url

    def sanitize(yml_str: str) -> str:
        """
        Scan the input string and sanitize any URLs that include
        credentials for source domains, returning the updated string.
        """
        # Match all URLs (i.e., `http://` and `https://` with or without credentials)
        url_regex = re.compile(r'https?://[^\s\'"]+')

        def replacer(match: re.Match) -> str:
            url = match.group(0)
            return sanitize_url(url)

        yml_sanitized_str = url_regex.sub(replacer, yml_str)

        return yml_sanitized_str

    env_manager = get_environment_manager()
    environment_cmd = get_environment_cmd()

    # Handle pixi separately since it writes a `pixi.lock` file
    if env_manager == "pixi":
        try:
            subprocess.run(
                environment_cmd,
                shell=True,
                check=True,
                stderr=subprocess.DEVNULL,
            )
            # https://pixi.sh/dev/reference/environment_variables/#environment-variables-set-by-pixi
            manifest_path = os.environ.get("PIXI_PROJECT_MANIFEST")
            lock_path = os.path.join(
                os.path.dirname(manifest_path) if manifest_path else os.getcwd(),
                "pixi.lock",
            )
            with open(lock_path, encoding="utf-8") as f:
                return sanitize(f.read())
        except Exception:
            return ""

    # For uv/conda/mamba environment managers, run the export command and process the output
    try:
        result = subprocess.run(
            environment_cmd,
            shell=True,
            check=True,
            stderr=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError:
        return ""

    raw_yml = result.stdout.strip()
    if not raw_yml:
        return ""

    if env_manager == "uv":
        return remove_comments(raw_yml) # Not sanitized, since uv doesn't use conda channels
    elif env_manager in ("conda", "mamba"):
        return sanitize(remove_metadata(raw_yml))

    raise RuntimeError(f"Unsupported environment manager: '{env_manager}'")


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


@parse_input_file_to_protocols_str.register(PackedPose)
@parse_input_file_to_protocols_str.register(Pose)
@parse_input_file_to_protocols_str.register(str)
def _parse_str(obj: str) -> List[str]:
    scores_dict = get_scores_dict(obj)
    return scores_dict["metadata"]["protocols"]


@singledispatch
def parse_input_file_to_instance_kwargs(obj: Any) -> NoReturn:
    raise InputFileError(obj)


@parse_input_file_to_instance_kwargs.register(PackedPose)
@parse_input_file_to_instance_kwargs.register(Pose)
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


def parse_init_file(
    input_file: str,
    input_packed_pose: Optional[PackedPose],
    skip_corrections: bool,
    init_from_file_kwargs: Dict[str, Any],
) -> Union[Tuple[Optional[PackedPose], PackedPose], NoReturn]:
    """
    Return a `tuple` object of the input `PackedPose` object and the output `PackedPose`
    object from a '.init' or '.init.bz2' file, verifying PyRosettaCluster metadata in the
    '.init' or '.init.bz2' file.
    """

    if not was_init_called():
        if skip_corrections != init_from_file_kwargs["skip_corrections"]:
            _skip_corrections_warning_msg = (
                "Please set the 'skip_corrections' keyword argument in the `reproduce()` function and in "
                "the 'init_from_file_kwargs' keyword arguments to the same value to silence this warning."
            )
            if skip_corrections and not init_from_file_kwargs["skip_corrections"]:
                warnings.warn(
                    "Skipping ScoreFunction corrections for the PyRosettaCluster task but "
                    + f"not for the host node PyRosetta initialization from the '.init' file! {_skip_corrections_warning_msg}",
                    UserWarning,
                    stacklevel=3,
                )
            elif not skip_corrections and init_from_file_kwargs["skip_corrections"]:
                warnings.warn(
                    "Skipping ScoreFunction corrections for the host node PyRosetta initialization "
                    + f"from the '.init' file but not for the PyRosettaCluster task! {_skip_corrections_warning_msg}",
                    UserWarning,
                    stacklevel=3,
                )
        try:
            io.init_from_file(input_file, **init_from_file_kwargs)
        except BufferError as ex:
            raise BufferError(
                f"{ex}. Please set a larger 'max_decompressed_bytes' parameter in the 'init_from_file_kwargs' "
                + "keyword argument of the `reproduce()` function to initialize PyRosetta with the input "
                + f"PyRosetta initialization file: '{input_file}'"
            )
        except Exception as ex:
            raise Exception(
                f"{type(ex).__name__}: {ex}. Could not initialize PyRosetta from the input PyRosetta initialization "
                + f"file '{input_file}' using `pyrosetta.init_from_file()` keyword arguments: '{init_from_file_kwargs}'. "
                + "Please ensure `pyrosetta.distributed.io.init_from_file()` runs with the '.init' or '.init.bz2' file "
                + "separately before passing it into `reproduce()`, and update any necessary `pyrosetta.init_from_file()` "
                + "keyword arguments in the 'init_from_file_kwargs' keyword argument parameter of the `reproduce()` function. "
                + "The '.init' or '.init.bz2' file may also be passed to `reproduce()` after a separate PyRosetta initialization."
            )
    else:
        _skip_corrections_warning_msg = (
            "Please ensure that PyRosetta was initialized from the same PyRosetta initialization file using "
            + f"`pyrosetta.init_from_file(skip_corrections={skip_corrections})` before running the `reproduce()` "
            + "function. To silence this warning, please ensure that PyRosetta is not already initialized before "
            + f"running the `reproduce()` function with the input PyRosetta initialization file: '{input_file}'"
        )
        if skip_corrections:
            warnings.warn(
                "Skipping ScoreFunction corrections for the PyRosettaCluster task but PyRosetta is already "
                + "initialized on the host node (with or without skipped ScoreFunction corrections)! "
                + _skip_corrections_warning_msg,
                UserWarning,
                stacklevel=3,
            )
        else:
            warnings.warn(
                "Preserving ScoreFunction corrections for the PyRosettaCluster task but PyRosetta is already "
                + "initialized on the host node (with or without preserved ScoreFunction corrections)! "
                + _skip_corrections_warning_msg,
                UserWarning,
                stacklevel=3,
            )

    _input_packed_pose, _output_packed_pose = get_poses_from_init_file(input_file, verify=True)
    if _output_packed_pose is None:
        raise ValueError(
            f"The input '.init' file does not contain an output decoy from a PyRosettaCluster simulation: '{input_file}'. "
            + "To reproduce from a '.init' file, please ensure that `pyrosetta.distributed.cluster.export_init_file()` "
            + "was run on the original decoy output file, or that the `PyRosettaCluster(output_decoy_types=['.init'])` "
            + "PyRosetta initialization file output decoy type was enabled in the original PyRosettaCluster simulation."
        )

    input_packed_pose = parse_input_packed_pose(input_packed_pose)
    if input_packed_pose is not None and not identical_b64_poses(input_packed_pose, _input_packed_pose):
        _input_packed_pose_error_msg = (
            "the input `PackedPose` object from the original PyRosettaCluster simulation is not identical "
            "to the provided 'input_packed_pose' keyword argument parameter of the `reproduce()` function"
        ) if _input_packed_pose is not None else (
            "the '.init' file does not contain an input `PackedPose` object from the original PyRosettaCluster simulation"
        )
        raise TypeError(
            "Please set the 'input_packed_pose' keyword argument parameter to a `NoneType` object when "
            + f"reproducing from a '.init' file, because {_input_packed_pose_error_msg}."
        )

    return (_input_packed_pose, _output_packed_pose)


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


def identical_b64_poses(
    packed_pose_1: Union[Pose, PackedPose],
    packed_pose_2: Union[Pose, PackedPose],
) -> bool:
    """Test whether b64-encoded pickled Pose objects are identical."""
    return io.to_base64(packed_pose_1) == io.to_base64(packed_pose_2)


def is_bytes(obj: Any) -> bool:
    return isinstance(obj, bytes)


def is_packed(obj: Any) -> bool:
    return isinstance(obj, PackedPose)


def is_dict(obj: Any) -> bool:
    return isinstance(obj, dict)
