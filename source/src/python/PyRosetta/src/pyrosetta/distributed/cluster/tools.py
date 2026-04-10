# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

try:
    import toolz
    from distributed import Client
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.tools' requires the "
        + "third-party packages 'distributed' and 'toolz' as dependencies!\n"
        + "Please install these packages into your virtual environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/distributed/\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import inspect
import json
import logging
import os
import tempfile
import warnings

from functools import wraps
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.exceptions import PyRosettaIsNotInitializedError
from pyrosetta.rosetta.basic import was_init_called

from pyrosetta.distributed.cluster.converters import _parse_protocols
from pyrosetta.distributed.cluster.converter_tasks import (
    is_dict,
    is_empty,
    get_protocols_list_of_str,
    parse_client,
    parse_decoy_name,
    parse_init_file,
    parse_input_file_to_instance_kwargs,
    parse_input_file_to_instance_metadata_kwargs,
    parse_instance_kwargs,
    parse_scorefile,
    reserve_scores_in_results,
)
from pyrosetta.distributed.cluster.core import PyRosettaCluster
from pyrosetta.distributed.cluster.io import secure_read_pickle
from pyrosetta.distributed.cluster.serialization import update_scores
from pyrosetta.distributed.cluster.type_defs import (
    Any,
    Callable,
    Dict,
    FloatOrInt,
    Generator,
    List,
    ListOrTuple,
    Optional,
    PoseOrPackedPose,
    PyRosettaProtocol,
    PyRosettaProtocolResults,
    PyRosettaProtocolType,    
    Tuple,
    Union,
    cast,
)


def get_protocols(
    protocols: Optional[Union[PyRosettaProtocol, ListOrTuple[PyRosettaProtocol]]]= None,
    input_file: Optional[Union[str, PoseOrPackedPose]] = None,
    scorefile: Optional[str] = None,
    decoy_name: Optional[str] = None,
) -> List[PyRosettaProtocol]:
    """
    Given an input file that was written by `PyRosettaCluster`, or a scorefile with full simulation records
    that was written by `PyRosettaCluster` and a decoy name, if an iterable of PyRosetta protocols is provided
    then validate the PyRosetta protocols against those in the full simulation record, otherwise if PyRosetta
    protocols are not provided then attempt to return the PyRosetta protocols from back two frames of the
    current scope where callable names may match the PyRosetta protocol names in the full simulation record.
    If the `protocols` keyword argument value is `None`, it is recommended to validate that the returned
    PyRosetta protocols are indeed those executed in the original `PyRosettaCluster` simulation.

    *Warning*: This function uses the `pickle` module to deserialize pickled `Pose` objects and pickled
    `pandas.DataFrame` objects. Using the `pickle` module is not secure, so please only run with input files
    you trust. Learn more about the `pickle` module and its security
    `here <https://docs.python.org/3/library/pickle.html>`_.

    Args:
        `protocols`: `PyRosettaProtocol | list[PyRosettaProtocol] | tuple[PyRosettaProtocol, ...] | None`
            An ordered iterable of callables (each of `types.GeneratorType` and/or `types.FunctionType` type) or
            a single callable of type `types.GeneratorType` or `types.FunctionType`, specifying the user-defined
            PyRosetta protocol(s) to execute for the reproduction simulation. If `None`, then the PyRosetta
            protocols are automatically detected (by a best effort only) from back two frames of the current
            scope.

            Default: `None`

        `input_file`: `str | Pose | PackedPose | None`
            A `str` object specifying the path to the ".pdb", ".pdb.bz2", ".pkl_pose", ".pkl_pose.bz2",
            ".b64_pose", ".b64_pose.bz2", ".init" or ".init.bz2" file from which to extract `PyRosettaCluster`
            instance attributes. If `input_file` is provided, then ignore the `scorefile` and `decoy_name`
            keyword arguments. Note that ".pkl_pose", ".pkl_pose.bz2", ".b64_pose", ".b64_pose.bz2", ".init"
            and ".init.bz2" files contain pickled `Pose` objects that are deserialized using the
            `SecureSerializerBase` class in PyRosetta upon calling the this function, but please still
            only input these file types if you know and trust their source. Learn more
            `here <https://docs.python.org/3/library/pickle.html>`_.

            Default: `None`

        `scorefile`: `str | None`
            A `str` object specifying the path to a JSON Lines (JSONL)-formatted scorefile or pickled
            `pandas.DataFrame` scorefile from a `PyRosettaCluster` simulation from which to extract
            `PyRosettaCluster` instance attributes. If `scorefile` is provided, then `decoy_name` must also be
            provided. In order to use a scorefile, it must contain full simulation records from the original
            `PyRosettaCluster` simulation; i.e., the `simulation_records_in_scorefile` keyword argument value
            was set to `True`. Note that in order to securely load pickled `pandas.DataFrame` objects, please
            ensure that `pyrosetta.secure_unpickle.add_secure_package("pandas")` has been run. If using `pandas`
            version `>=3.0.0`, PyArrow-backed datatypes may be enabled by default; in this case, please ensure
            that `pyrosetta.secure_unpickle.add_secure_package("pyarrow")` has also first been run.

            Default: `None`

        `decoy_name`: `str | None`
            A `str` object specifying the decoy name for which to extract `PyRosettaCluster` instance
            attributes. If `decoy_name` is provided, then `scorefile` must also be provided.

            Default: `None`

    Returns:
        A `list` object of callable PyRosetta protocols as determined from the full simulation record. If the
        `protocols` keyword argument value is `None`, then attempt to return the callable PyRosetta protocols
        (by a best effort only) from back two frames of the current scope where callable names may match the
        PyRosetta protocol names in the full simulation record.
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
            "The original user-defined PyRosetta protocols list and the `protocols` keyword argument "
            + "value have different lengths! Cannot reproduce!"
        )
        if original_protocols_list_of_str == input_protocols_list_of_str:
            logging.info(
                "The `protocols` keyword argument value matches the user-defined PyRosetta protocols "
                + "from the original production simulation. Continuing with the reproduction simulation."
            )
        for i in range(len(original_protocols_list_of_str)):
            if original_protocols_list_of_str[i] != input_protocols_list_of_str[i]:
                logging.warning(
                    f"The original user-defined PyRosetta protocol '{original_protocols_list_of_str[i]}' "
                    + f"appears to have changed names to '{input_protocols_list_of_str[i]}'! "
                    + f"Please verify that '{input_protocols_list_of_str[i]}' is functionally equivalent to "
                    + f"'{original_protocols_list_of_str[i]}', otherwise the reproduction simulation "
                    + "will not reproduce the original decoy! Continuing with the reproduction simulation."
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
                    f"Automatically detected user-defined PyRosetta protocol '{protocol_name}' "
                    + "in the local variables of the current frame."
                )
                protocols.append(scope.f_locals[protocol_name])
            elif protocol_name in scope.f_globals:
                logging.info(
                    f"Automatically detected user-defined PyRosetta protocol '{protocol_name}' "
                    + "in the global variables of the current frame."
                )
                protocols.append(scope.f_globals[protocol_name])
            else:
                raise RuntimeError(
                    f"The original user-defined PyRosetta protocol '{protocol_name}' "
                    + "could not be found back two frames of the current scope! Please verify that the "
                    + "original PyRosetta protocols are defined in the same scope as they were in "
                    + "the original production run. Alternatively, you may pass the original "
                    + "PyRosetta protocols (defined under the new scope) maintaining "
                    + "the original protocol order as a value to the `protocols` keyword argument."
                )

    return _parse_protocols(protocols)


def get_instance_kwargs(
    input_file: Optional[Union[str, PoseOrPackedPose]] = None,
    scorefile: Optional[str] = None,
    decoy_name: Optional[str] = None,
    skip_corrections: Optional[bool] = None,
    with_metadata_kwargs: Optional[bool] = None,
) -> Union[
    Dict[str, Any],
    Tuple[Dict[str, Any], Dict[str, Any]]
]:
    """
    Given an input file that was written by `PyRosettaCluster`, or a scorefile with full simulation records
    that was written by `PyRosettaCluster` and a decoy name, return the `PyRosettaCluster` instance attributes
    (and optionally the "metadata" keyword arguments) needed to reproduce the decoy using `PyRosettaCluster`.

    *Warning*: This function uses the `pickle` module to deserialize pickled `Pose` objects and pickled
    `pandas.DataFrame` objects. Using the `pickle` module is not secure, so please only run with input files
    you trust. Learn more about the `pickle` module and its security
    `here <https://docs.python.org/3/library/pickle.html>`_.

    Args:
        `input_file`: `str | Pose | PackedPose | None`
            A `str` object specifying the path to the ".pdb", ".pdb.bz2", ".pkl_pose", ".pkl_pose.bz2",
            ".b64_pose", ".b64_pose.bz2", ".init" or ".init.bz2" file from which to extract `PyRosettaCluster`
            instance attributes. If `input_file` is provided, then ignore the `scorefile` and `decoy_name`
            keyword arguments. Note that ".pkl_pose", ".pkl_pose.bz2", ".b64_pose", ".b64_pose.bz2", ".init"
            and ".init.bz2" files contain pickled `Pose` objects that are deserialized using the
            `SecureSerializerBase` class in PyRosetta upon calling the this function, but please still
            only input these file types if you know and trust their source. Learn more
            `here <https://docs.python.org/3/library/pickle.html>`_.

            Default: `None`

        `scorefile`: `str | None`
            A `str` object specifying the path to a JSON Lines (JSONL)-formatted scorefile or pickled
            `pandas.DataFrame` scorefile from a `PyRosettaCluster` simulation from which to extract
            `PyRosettaCluster` instance attributes. If `scorefile` is provided, then `decoy_name` must also be
            provided. In order to use a scorefile, it must contain full simulation records from the original
            `PyRosettaCluster` simulation; i.e., the `simulation_records_in_scorefile` keyword argument value
            was set to `True`. Note that in order to securely load pickled `pandas.DataFrame` objects, please
            ensure that `pyrosetta.secure_unpickle.add_secure_package("pandas")` has been run. If using `pandas`
            version `>=3.0.0`, PyArrow-backed datatypes may be enabled by default; in this case, please ensure
            that `pyrosetta.secure_unpickle.add_secure_package("pyarrow")` has also first been run.

            Default: `None`

        `decoy_name`: `str | None`
            A `str` object specifying the decoy name for which to extract `PyRosettaCluster` instance
            attributes. If `decoy_name` is provided, then `scorefile` must also be provided.

            Default: `None`

        `skip_corrections`: `bool | None`
            A `bool` object specifying whether or not to skip any `ScoreFunction` corrections specified in the
            `PyRosettaCluster` task's PyRosetta initialization options (extracted from the full simulation
            record in the `input_file` or `scorefile` keyword argument value). If `None`, then `False`.

            Default: `None`

        `with_metadata_kwargs`: `bool | None`
            A `bool` object specifying whether or not to return a `tuple` object with the `PyRosettaCluster`
            instance attributes as the first element and the "metadata" keyword arguments as the second element.
            If `None`, then `False`.

            Default: `None`

    Returns:
        A `dict` object of `PyRosettaCluster` instance attributes, or a `tuple` object of `dict` objects with
        the `PyRosettaCluster` instance attributes as the first element and the "metadata" keyword arguments as
        the second element when the `with_metadata_kwargs` keyword argument value is set to `True`.
    """

    _simulation_records_in_scorefile_msg = (
        "The `scorefile` keyword argument value does not contain the full simulation records. "
        + "In order to reproduce a decoy using a scorefile, the `PyRosettaCluster` instance attribute "
        + "'simulation_records_in_scorefile' must have been set to `True` in the original simulation. "
        + "Please provide an output decoy file that was written by `PyRosettaCluster` to the `input_file` "
        + "keyword argument value, or an output scorefile with full simulation records that was written by "
        + "`PyRosettaCluster` to the `scorefile` keyword argument value, in order to reproduce."
    )
    if input_file:
        if scorefile or decoy_name:
            warnings.warn(
                "Received `input_file` and either `scorefile` or `decoy_name` keyword arguments. "
                + "Ignoring `scorefile` and `decoy_name` values and using the `input_file` value!",
                UserWarning,
                stacklevel=2,
            )
        if with_metadata_kwargs:
            instance_kwargs, metadata_kwargs = parse_input_file_to_instance_metadata_kwargs(input_file)
        else:
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
                            "Received a `scorefile` keyword argument value that does not appear "
                            + f"to be JSONL-formatted: '{scorefile}'"
                        )
                    if all(k in scorefile_entry for k in ("metadata", "instance")):
                        if "decoy_name" in scorefile_entry["metadata"]:
                            if scorefile_entry["metadata"]["decoy_name"] == decoy_name:
                                instance_kwargs = scorefile_entry["instance"]
                                if with_metadata_kwargs:
                                    metadata_kwargs = scorefile_entry["metadata"]
                                break
                    else:
                        raise NotImplementedError(_simulation_records_in_scorefile_msg)
        else:
            try:
                df = secure_read_pickle(scorefile, compression="infer")
            except:
                raise TypeError(
                    "Received a `scorefile` keyword argument value that does not appear to be readable by "
                    + f"`pyrosetta.distributed.cluster.io.secure_read_pickle(compression='infer')`: '{scorefile}'"
                )
            if all(k in df.columns for k in ("metadata", "instance")):
                for instance, metadata in df[["instance", "metadata"]].values:
                    if "decoy_name" in metadata:
                        if metadata["decoy_name"] == decoy_name:
                            instance_kwargs = dict(instance)
                            if with_metadata_kwargs:
                                metadata_kwargs = dict(metadata)
                            break
            else:
                raise NotImplementedError(_simulation_records_in_scorefile_msg)
        if instance_kwargs is None:
            raise KeyError(
                f"The `decoy_name` keyword argument value is not in the provided scorefile: '{scorefile}'"
            )
    else:
        raise NotImplementedError(
            "Either an `input_file` keyword argument value, or both `scorefile` and `decoy_name` keyword "
            + "argument values, must be provided."
        )
    assert isinstance(
        instance_kwargs, dict
    ), f"Returned instance keyword arguments are not of type `dict`: {type(instance_kwargs)}"
    if with_metadata_kwargs:
        assert isinstance(
            metadata_kwargs, dict
        ), f"Returned metadata keyword arguments are not of type `dict`: {type(metadata_kwargs)}"

    if skip_corrections:
        assert isinstance(
            instance_kwargs["tasks"], dict
        ), f"The 'tasks' instance attribute must be of type `dict`: {type(instance_kwargs['tasks'])}"
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
                            "Cannot skip `ScoreFunction` corrections because the original `PyRosettaCluster` simulation "
                            + "did not output results with normalized Rosetta command-line options or configure "
                            + "the task's Rosetta command-line options as an instance of `dict`. Please disable the "
                            + "'skip_corrections' keyword argument and try again."
                        )
                else:
                    raise TypeError(
                        f"The `PyRosettaCluster` task key '{option}' must have a value of type `dict` or `str`. "
                        + f"Received: {type(instance_kwargs['tasks'][option])}"
                    )

    if with_metadata_kwargs:
        return instance_kwargs, metadata_kwargs
    else:
        return instance_kwargs


def reserve_scores(func: PyRosettaProtocolType) -> PyRosettaProtocolType:
    """
    A decorator for any user-defined PyRosetta protocol. If any non-scoreterm keys are present in the input
    `PackedPose` object's `Pose.cache` dictionary, and if they are deleted during execution of the decorated
    PyRosetta protocol, then restore those non-scoreterm keys and their values back into the `Pose.cache`
    dictionary after execution. If any keys are present in the input `PackedPose` object's `Pose.cache`
    dictionary and also present in the returned or yielded output `Pose` or `PackedPose` object's `Pose.cache`
    dictionary, then do not set the original values back into the `Pose.cache` dictionary after execution (that
    is, do not overwrite the new values in the `Pose.cache` dictionary). Any new keys acquired in the decorated
    PyRosetta protocol will never be overwritten. This allows users to maintain keys and values acquired in
    upstream PyRosetta protocols if needing to execute RosettaScripts `Mover` objects that happen to delete
    cached scores from `Pose` objects. Note that this decorator reserves keys and values from the input
    `PackedPose.scores` and `Pose.cache` dictionaries, and any keys and values restored to any output
    `Pose.cache` dictionaries are set as `SimpleMetrics` metrics.

    *Warning*: This decorator uses the `pickle` module to deserialize pickled `Pose` objects and arbitrary
    Python types in `Pose.cache` dictionary. Using the `pickle` module is not secure, so please only run with
    input data you trust. Learn more about the `pickle` module and its security
    `here <https://docs.python.org/3/library/pickle.html>`_.

    Example:

        >>> @reserve_scores
        ... def my_pyrosetta_protocol(packed_pose, /, **kwargs):
        ...     from pyrosetta import MyMover
        ...     pose = packed_pose.pose
        ...     MyMover().apply(pose) # Deletes scores
        ...     return pose

    Args:
        `func`: `PyRosettaProtocol`
            A callable of type `types.GeneratorType` or `types.FunctionType` representing a user-defined
            PyRosetta protocol.

    Returns:
        The output from the decorated user-defined PyRosetta protocol with reserved scores from the input
        `PackedPose` object's `Pose.cache` dictionary updated into the output `Pose` or `PackedPose`
        object's (or objects') `Pose.cache` dictionary (or dictionaries).
    """
    import pyrosetta  # noqa
    import pyrosetta.distributed  # noqa

    @wraps(func)
    def wrapper(packed_pose: Optional[PackedPose], **kwargs: Any) -> PyRosettaProtocolResults:
        """Wrapper function for the `reserve_scores` decorator."""

        if packed_pose is not None:
            _reserved_pose = update_scores(packed_pose).pose
        else:
            _reserved_pose = None
        _output = func(packed_pose, **kwargs)
        # Only deserialize after the user-provided PyRosetta protocol finished executing, giving
        # the user an opportunity to add secure packages to the unpickle-allowed list as necessary
        _scores_dict = dict(_reserved_pose.cache) if _reserved_pose is not None else {}

        return reserve_scores_in_results(_output, _scores_dict, func.__name__)

    return cast(PyRosettaProtocolType, wrapper)


def requires_packed_pose(func: PyRosettaProtocolType) -> PyRosettaProtocolType:
    """
    A decorator for any user-defined PyRosetta protocol. If a PyRosetta protocol requires that the first
    positional-or-keyword parameter be a non-empty `PackedPose` object, then immediately return any bound empty
    `PackedPose` or `None` objects and skip the decorated PyRosetta protocol, otherwise run the decorated
    PyRosetta protocol. If using `PyRosettaCluster(filter_results=False)` and the preceding PyRosetta protocol
    produces either `None`, an empty `Pose` object, or an empty `PackedPose` object, then an empty `PackedPose`
    object is distributed to the next PyRosetta protocol, in which case the next protocol and/or any downstream
    PyRosetta protocols are skipped if they are decorated with this decorator. If using
    `PyRosettaCluster(ignore_errors=True)` and a standard Python exception is raised or a Rosetta segmentation
    fault is thrown in the preceding PyRosetta protocol, then `None` is distributed to the next PyRosetta
    protocol, in which case the next PyRosetta protocol and/or any downstream PyRosetta protocols are skipped if
    they are decorated with this decorator.

    Example:

        >>> @requires_packed_pose
        ... def my_pyrosetta_protocol(packed_pose, /, **kwargs):
        ...     assert not packed_pose.empty() or packed_pose.pose.size() > 0
        ...     return packed_pose

    Args:
        `func`: `PyRosettaProtocol`
            A callable of type `types.GeneratorType` or `types.FunctionType` representing a user-defined
            PyRosetta protocol.

    Returns:
        The first positional-or-keyword parameter if it is an empty `PackedPose` object or `None`, otherwise
        the produced results from the decorated user-defined PyRosetta protocol.
    """

    @wraps(func)
    def wrapper(packed_pose: Optional[PackedPose], **kwargs: Any) -> PyRosettaProtocolResults:
        """Wrapper function for the `requires_packed_pose` decorator."""

        _msg = "User-provided PyRosetta protocol '{0}' received and is duly returning {1} object."
        if is_empty(packed_pose):
            logging.info(_msg.format(func.__name__, "an empty `PackedPose`"))
            return packed_pose
        elif packed_pose is None:
            logging.info(_msg.format(func.__name__, "a `NoneType`"))
            return packed_pose
        else:
            return func(packed_pose, **kwargs)

    return cast(PyRosettaProtocolType, wrapper)


def reproduce(
    input_file: Optional[str] = None,
    scorefile: Optional[str] = None,
    decoy_name: Optional[str] = None,
    protocols: Optional[Union[PyRosettaProtocol, ListOrTuple[PyRosettaProtocol]]] = None,
    client: Optional[Client] = None,
    clients: Optional[ListOrTuple[Client]] = None,
    input_packed_pose: Optional[PoseOrPackedPose] = None,
    instance_kwargs: Optional[Dict[str, Any]] = None,
    clients_indices: Optional[ListOrTuple[int]] = None,
    resources: Optional[ListOrTuple[Optional[Dict[str, FloatOrInt]]]] = None,
    retries: Optional[Union[int, ListOrTuple[int]]] = None,
    skip_corrections: bool = False,
    init_from_file_kwargs: Optional[Dict[str, Any]] = None,
) -> None:
    """
    Given an input file (or a scorefile with full simulation records and a decoy name) that was written by
    `PyRosettaCluster` and any additional `PyRosettaCluster` class keyword arguments, execute the decoy
    reproduction simulation with a new instance of `PyRosettaCluster`.

    *Warning*: This method uses the `cloudpickle` and `pickle` modules to serialize and deserialize `Pose`
    objects, arbitrary Python types in `Pose.cache` dictionaries, `pandas.DataFrame` objects (if configured),
    user-defined task dictionaries, user-defined PyRosetta protocols, and other user-provided data. Using the
    `cloudpickle` and `pickle` modules is not secure, so please only run this method with input data you fully
    understand and trust. Learn more about the `cloudpickle` and `pickle` modules and their security
    `here <https://github.com/cloudpipe/cloudpickle>`_ and
    `here <https://docs.python.org/3/library/pickle.html>`_.

    Args:
        `input_file`: `str | None`
            A `str` object specifying the path to the ".pdb", ".pdb.bz2", ".pkl_pose", ".pkl_pose.bz2",
            ".b64_pose", ".b64_pose.bz2", ".init" or ".init.bz2" file from which to extract `PyRosettaCluster`
            instance attributes. If `input_file` is provided, then ignore the `scorefile` and `decoy_name`
            keyword arguments. If a ".init" or ".init.bz2" file is provided and PyRosetta is not yet
            initialized, then first initialize PyRosetta with the PyRosetta initialization file (see the
            `init_from_file_kwargs` keyword argument). Note that ".pkl_pose", ".pkl_pose.bz2", ".b64_pose",
            ".b64_pose.bz2", ".init" and ".init.bz2" files contain pickled `Pose` objects that are deserialized
            using the `SecureSerializerBase` class in PyRosetta upon calling the `reproduce` function, but
            please still only input these file types if you know and trust their source. Learn more
            `here <https://docs.python.org/3/library/pickle.html>`_.

            Default: `None`

        `scorefile`: `str | None`
            A `str` object specifying the path to a JSON Lines (JSONL)-formatted scorefile or pickled
            `pandas.DataFrame` scorefile from a `PyRosettaCluster` simulation from which to extract
            `PyRosettaCluster` instance attributes. If `scorefile` is provided, then `decoy_name` must also be
            provided. In order to use a scorefile, it must contain full simulation records from the original
            `PyRosettaCluster` simulation; i.e., the `simulation_records_in_scorefile` keyword argument value
            was set to `True`. Note that in order to securely load pickled `pandas.DataFrame` objects, please
            ensure that `pyrosetta.secure_unpickle.add_secure_package("pandas")` has been run. If using `pandas`
            version `>=3.0.0`, PyArrow-backed datatypes may be enabled by default; in this case, please ensure
            that `pyrosetta.secure_unpickle.add_secure_package("pyarrow")` has also first been run.

            Default: `None`

        `decoy_name`: `str | None`
            A `str` object specifying the decoy name for which to extract `PyRosettaCluster` instance
            attributes. If `decoy_name` is provided, then `scorefile` must also be provided.

            Default: `None`

        `protocols`: `PyRosettaProtocol | list[PyRosettaProtocol] | tuple[PyRosettaProtocol, ...] | None`
            An ordered iterable of callables (each of `types.GeneratorType` and/or `types.FunctionType` type) or
            a single callable of type `types.GeneratorType` or `types.FunctionType`, specifying the user-defined
            PyRosetta protocol(s) to execute for the reproduction simulation. If `None`, the PyRosetta protocols
            are automatically detected (by a best effort only) in the scope from which the `reproduce` function
            is called.

            Default: `None`

        `client`: `distributed.Client | None`
            An initialized Dask `distributed.Client` object to be used as the Dask client interface to the local
            or remote Dask cluster. If `None`, then `PyRosettaCluster` initializes its own Dask client based on
            the `scheduler` keyword argument value (see `instance_kwargs` keyword argument). Deprecated by the
            `clients` keyword argument, but supported for legacy purposes. Either or both of the `client` or
            `clients` keyword argument values must be `None`.

            Default: `None`

        `clients`: `list[distributed.Client] | tuple[distributed.Client, ...] | None`
            A `list` or `tuple` object of initialized Dask `distributed.Client` objects to be used as the Dask
            client interface(s) to the local or remote Dask cluster(s). If `None`, then `PyRosettaCluster`
            initializes its own Dask client based on the `scheduler` keyword argument value (see
            `instance_kwargs` keyword argument). Optionally used in combination with the `clients_indices`
            keyword argument. Either or both of the `client` or `clients` keyword argument values must be
            `None`.

            Default: `None`

        `input_packed_pose`: `Pose | PackedPose | None`
            An input `PackedPose` object that is accessible via the first positional-or-keyword parameter of the
            first user-defined PyRosetta protocol.

            Default: `None`

        `instance_kwargs`: `dict[str, Any] | None`
            A `dict` object of valid `PyRosettaCluster` keyword arguments which will override any
            `PyRosettaCluster` instance attributes that were used to generate the original decoy and that were
            stored in the full simulation record.

            Default: `None`

        `clients_indices`: `list[int] | tuple[int, ...] | None`
            A `list` or `tuple` object of `int` objects, where each `int` object represents a zero-based index
            corresponding to the initialized Dask `distributed.Client` object(s) passed to the `clients`
            keyword argument value. If not `None`, then the length of the `clients_indices` object must equal
            the number of PyRosetta protocols (see `protocols` keyword argument).

            Default: `None`

        `resources`: `list[dict[str, float | int] | None] | tuple[dict[str, float | int] | None, ...] | None`
            A `list` or `tuple` object of `dict` objects, where each `dict` object represents an abstract,
            arbitrary resource to constrain which Dask workers execute the user-defined PyRosetta protocols. If
            `None`, then do not impose resource constaints on any PyRosetta protocols. If not `None`, then the
            length of the `resources` object must equal the number of PyRosetta protocols passed to the
            `protocols` keyword argument, such that each resource specified indicates the unique resource
            constraints for the protocol at the corresponding index of the PyRosetta protocols passed to the
            `protocols` keyword argument. Note that this feature is only useful when one passes in their own
            instantiated Dask client(s) with Dask workers set up with various resource constraints. If Dask
            workers were not instantiated to satisfy the specified resource constraints, PyRosetta protocols
            will hang indefinitely by design because the Dask scheduler is waiting for Dask workers that meet
            the specified resource constraints so that it may schedule these tasks. Unless Dask workers were
            created with these resource tags applied, the PyRosetta protocols will not run.

            See https://distributed.dask.org/en/stable/resources.html for more information.

            Default: `None`

        `retries`: `list[int] | tuple[int, ...] | int | None`
            A `list` or `tuple` of `int` objects, where each `int` object (≥0) sets the number of allowed
            automatic retries of each failed task that was applied to the corresponding user-defined PyRosetta
            protocol (i.e., indexed the same as `client_indices` keyword argument value). If an `int` object
            (≥0) is provided, then apply that number of allowed automatic retries to all PyRosetta protocols.
            If `None` is provided, then no explicit retries are allowed. If not `None` and not an `int` object,
            then the length of this value must equal the number of PyRosetta protocols passed to the
            `PyRosettaCluster.distribute` method, and each `int` value determines the number of automatic
            retries the Dask scheduler allows for that the tasks applied to that PyRosetta protocol. Allowing
            retries of failed tasks may be useful if the PyRosetta protocol raises a standard Python exception
            or Rosetta throws a segmentation fault in the `billiard` subprocess while the Dask worker remains
            alive and the value of the `ignore_errors` key in the `instance_kwargs` keyword argument value is
            set to `False`. If `ignore_errors` is set to `True`, then protocols failing due to standard Python
            exceptions or Rosetta segmentation faults will still be considered successes, and this keyword
            argument has no effect since these PyRosetta protocol errors are ignored. Note that if a compute
            resource executing a PyRosetta protocol is preempted, then the Dask worker process does not remain
            alive and the Dask scheduler registers that failed task as incomplete or cancelled. In this case,
            the number of allowed task retries is controlled by the Dask configuration parameter
            `distributed.scheduler.allowed-failures`; please use the `max_task_replices` and `task_registry`
            keys of the dictionary passed to the `instance_kwargs` keyword argument value for further
            configuration of task retries after compute resource preemption.

            See https://distributed.dask.org/en/latest/scheduling-state.html#task-state for more information.

            Default: `None`

        `skip_corrections`: `bool`
            A `bool` object specifying whether or not to skip any `ScoreFunction` corrections specified in the
            values of the `"options"` or `"extra_options"` keys from task dictionary of the original simulation.
            (extracted from the full simulation record from either the `input_file` or `scorefile` keyword
            argument value), which are set in-code upon PyRosetta initialization. If the current PyRosetta build
            and Conda/Mamba/uv/Pixi environment are identical to those used for the original simulation, this
            keyword argument value may be set to `True` to enable the reproduced output decoy file to be used
            for successive reproductions; however, some `ScoreFunction` corrections may be critical for decoy
            reproducibility. Therefore, it may be prudent to set this keyword argument value to `False` unless
            there is absolute certainty that `ScoreFunction` corrections may be skipped without influencing
            decoy reproducibility. If reproducing from a PyRosetta initialization file, it is recommended to also
            set the value of the `skip_corrections` key from the `init_from_file_kwargs` keyword argument value
            to the same boolean value.

            See https://docs.rosettacommons.org/docs/latest/full-options-list#corrections and
            https://docs.rosettacommons.org/docs/latest/full-options-list#1-corrections for more information.

            Default: `False`

        `init_from_file_kwargs`: `dict[str, Any] | None`
            A `dict` object to override the default `pyrosetta.init_from_file` keyword arguments if the
            `input_file` keyword argument value is a path to a PyRosetta initialization file, otherwise it is
            has no effect. See the
            `pyrosetta.init_from_file <pyrosetta.utility.initialization.PyRosettaInitFileParser.init_from_file>`_
            docstring for more information.

            Default: A `dict` object with the following entries:
                `output_dir`:
                    Default: `os.path.join(tempfile.TemporaryDirectory(prefix="PyRosettaCluster_reproduce_").name, "pyrosetta_init_input_files")`
                `skip_corrections`:
                    Default: The `skip_corrections` keyword argument value from the `reproduce` function.
                `relative_paths`:
                    Default: `True`
                `dry_run`:
                    Default: `False`
                `max_decompressed_bytes`:
                    Default: `pow(2, 30)` or 1 GiB.
                `restore_rg_state`:
                    Default: `True`
                `database`:
                    Default: `None`
                `verbose`:
                    Default: `True`
                `set_logging_handler`:
                    Default: `"logging"`
                `notebook`:
                    Default: `None`
                `silent`:
                    Default: `False`

    Returns:
        `None`
    """

    if not isinstance(skip_corrections, bool):
        raise TypeError(
            "The `skip_corrections` keyword argument value must be of type `bool`. "
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
                restore_rg_state=True,
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
                    + "keyword argument, please ensure `pyrosetta.init` or `pyrosetta.init_from_file` has been "
                    + "properly called (with the same residue type set as that used to generate the original '.pkl_pose', "
                    + "'.pkl_pose.bz2', '.b64_pose', or '.b64_pose.bz2' file) before running the `reproduce` function. "
                    + "If an output '.init' file from the original simulation is available, it is recommended to run "
                    + "`pyrosetta.init_from_file` with that '.init' file before running the `reproduce` function."
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
                    with_metadata_kwargs=False,
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
        priorities=None,
        retries=retries,
    )
    if isinstance(_tmp_dir, tempfile.TemporaryDirectory):
        _tmp_dir.cleanup()


def produce(**kwargs: Any) -> None:
    """
    A `PyRosettaCluster.distribute` method shim requiring the `protocols` keyword argument, and optionally any
    `PyRosettaCluster` keyword arguments or `PyRosettaCluster.distribute` keyword arguments.

    *Warning*: This method uses the `cloudpickle` and `pickle` modules to serialize and deserialize `Pose`
    objects, arbitrary Python types in `Pose.cache` dictionaries, `pandas.DataFrame` objects (if configured),
    user-defined task dictionaries, user-defined PyRosetta protocols, and other user-provided data. Using the
    `cloudpickle` and `pickle` modules is not secure, so please only run this method with input data you fully
    understand and trust. Learn more about the `cloudpickle` and `pickle` modules and their security
    `here <https://github.com/cloudpipe/cloudpickle>`_ and
    `here <https://docs.python.org/3/library/pickle.html>`_.

    Args:
        `**kwargs`: `Any`
            See the `PyRosettaCluster` docstring. The keyword arguments must also include `protocols`, an
            iterable object of callable or generator function objects specifying an ordered sequence of
            user-defined PyRosetta protocols to execute for the simulation (see `PyRosettaCluster.distribute`
            docstring). The keyword arguments may also optionally include `clients_indices`, `resources`,
            `priorities`, and `retries` (see `PyRosettaCluster.distribute` docstring).
    """
    # See `produce.__doc__` updated below

    protocols = kwargs.pop("protocols", None)
    clients_indices = kwargs.pop("clients_indices", None)
    resources = kwargs.pop("resources", None)
    priorities = kwargs.pop("priorities", None)
    retries = kwargs.pop("retries", None)
    PyRosettaCluster(**kwargs).distribute(
        protocols=protocols,
        clients_indices=clients_indices,
        resources=resources,
        priorities=priorities,
        retries=retries,
    )


run: Callable[..., None] = produce


@wraps(produce, assigned=("__doc__",), updated=())
def iterate(**kwargs: Any) -> Generator[Tuple[Optional[PackedPose], Dict[str, Any]], None, None]:
    # See assigned `iterate.__doc__` updated below
    protocols = kwargs.pop("protocols", None)
    clients_indices = kwargs.pop("clients_indices", None)
    resources = kwargs.pop("resources", None)
    priorities = kwargs.pop("priorities", None)
    retries = kwargs.pop("retries", None)
    for result in PyRosettaCluster(**kwargs).generate(
        protocols=protocols,
        clients_indices=clients_indices,
        resources=resources,
        priorities=priorities,
        retries=retries,
    ):
        yield result


produce.__doc__ += """
    Returns:
        `None`
    """


iterate.__doc__ = iterate.__doc__.replace(
    "PyRosettaCluster.distribute", "PyRosettaCluster.generate"
) + """
    Yields:
        ``(`PackedPose`, `dict`)`` tuples from the most recently executed user-defined PyRosetta protocol
        if `PyRosettaCluster(save_all=True)` is used, otherwise from the final PyRosetta protocol.
    """
