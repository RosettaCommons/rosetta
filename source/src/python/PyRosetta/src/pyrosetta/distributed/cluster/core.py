"""
PyRosettaCluster is a class for reproducible, high-throughput job distribution
of user-defined PyRosetta protocols efficiently parallelized on the user's
local computer, high-performance computing (HPC) cluster, or elastic cloud
computing infrastructure with available compute resources.

Args:
    tasks: A `list` of `dict` objects, a callable or called function returning
        a `list` of `dict` objects, or a callable or called generator yielding
        a `list` of `dict` objects. Each dictionary object element of the list
        is accessible via kwargs in the user-defined PyRosetta protocols.
        In order to initialize PyRosetta with user-defined PyRosetta command line
        options at the start of each user-defined PyRosetta protocol, either
        `extra_options` and/or `options` must be a key of each dictionary object,
        where the value is a `str`, `tuple`, `list`, `set`, or `dict` of
        PyRosetta command line options.
        Default: [{}]
    input_packed_pose: Optional input `PackedPose` object that is accessible via
        the first argument of the first user-defined PyRosetta protocol.
        Default: None
    seeds: A `list` of `int` objects specifying the random number generator seeds
        to use for each user-defined PyRosetta protocol. The number of seeds
        provided must be equal to the number of user-defined input PyRosetta
        protocols. Seeds are used in the same order that the user-defined PyRosetta
        protocols are executed.
        Default: None
    decoy_ids: A `list` of `int` objects specifying the decoy numbers to keep after
        executing user-defined PyRosetta protocols. User-provided PyRosetta
        protocols may return a list of `Pose` and/or `PackedPose` objects, or
        yield multiple `Pose` and/or `PackedPose` objects. To reproduce a
        particular decoy generated via the chain of user-provided PyRosetta
        protocols, the decoy number to keep for each protocol may be specified,
        where other decoys are discarded. Decoy numbers use zero-based indexing,
        so `0` is the first decoy generated from a particular PyRosetta protocol.
        The number of decoy_ids provided must be equal to the number of
        user-defined input PyRosetta protocols, so that one decoy is saved for each
        user-defined PyRosetta protocol. Decoy ids are applied in the same order
        that the user-defined PyRosetta protocols are executed.
        Default: None
    client: An initialized dask `distributed.client.Client` object to be used as
        the dask client interface to the local or remote compute cluster. If `None`,
        then PyRosettaCluster initializes its own dask client based on the
        `PyRosettaCluster(scheduler=...)` class attribute. Deprecated by the
        `PyRosettaCluster(clients=...)` class attribute, but supported for legacy
        purposes. Either or both of the `client` or `clients` attribute parameters
        must be `None`.
        Default: None
    clients: A `list` or `tuple` object of initialized dask `distributed.client.Client`
        objects to be used as the dask client interface(s) to the local or remote compute
        cluster(s). If `None`, then PyRosettaCluster initializes its own dask client based
        on the `PyRosettaCluster(scheduler=...)` class attribute. Optionally used in
        combination with the `PyRosettaCluster().distribute(clients_indices=...)` method.
        Either or both of the `client` or `clients` attribute parameters must be `None`.
        See the `PyRosettaCluster().distribute()` method docstring for usage examples.
        Default: None
    scheduler: A `str` of either "sge" or "slurm", or `None`. If "sge", then
        PyRosettaCluster schedules jobs using `SGECluster` with `dask-jobqueue`.
        If "slurm", then PyRosettaCluster schedules jobs using `SLURMCluster` with
        `dask-jobqueue`. If `None`, then PyRosettaCluster schedules jobs using
        `LocalCluster` with `dask.distributed`. If `PyRosettaCluster(client=...)`
        or `PyRosettaCluster(clients=...)` is provided, then 
        `PyRosettaCluster(scheduler=...)` is ignored.
        Default: None
    cores: An `int` object specifying the total number of cores per job, which
        is input to the `dask_jobqueue.SLURMCluster(cores=...)` argument or
        the `dask_jobqueue.SGECluster(cores=...)` argument.
        Default: 1
    processes: An `int` object specifying the total number of processes per job,
        which is input to the `dask_jobqueue.SLURMCluster(processes=...)` argument
        or the `dask_jobqueue.SGECluster(processes=...)` argument.
        This cuts the job up into this many processes.
        Default: 1
    memory: A `str` object specifying the total amount of memory per job, which
        is input to the `dask_jobqueue.SLURMCluster(memory=...)` argument or
        the `dask_jobqueue.SGECluster(memory=...)` argument.
        Default: "4g"
    scratch_dir: A `str` object specifying the path to a scratch directory where
        dask litter may go.
        Default: "/temp" if it exists, otherwise the current working directory
    min_workers: An `int` object specifying the minimum number of workers to
        which to adapt during parallelization of user-provided PyRosetta protocols.
        Default: 1
    max_workers: An `int` object specifying the maximum number of workers to
        which to adapt during parallelization of user-provided PyRosetta protocols.
        Default: 1000 if the initial number of `tasks` is <1000, else use the
            the initial number of `tasks`
    dashboard_address: A `str` object specifying the port over which the dask
        dashboard is forwarded. Particularly useful for diagnosing PyRosettaCluster
        performance in real-time.
        Default=":8787"
    nstruct: An `int` object specifying the number of repeats of the first
        user-provided PyRosetta protocol. The user can control the number of
        repeats of subsequent user-provided PyRosetta protocols via returning
        multiple clones of the output pose(s) from a user-provided PyRosetta
        protocol run earlier, or cloning the input pose(s) multiple times in a
        user-provided PyRosetta protocol run later.
        Default: 1
    compressed: A `bool` object specifying whether or not to compress the output
        ".pdb", ".pkl_pose", ".b64_pose", and ".init" files with `bzip2`, resulting
        in appending ".bz2" to decoy output files and PyRosetta initialization files.
        Also see the 'output_decoy_types' and 'output_init_file' keyword arguments.
        Default: True
    compression: A `str` object of 'xz', 'zlib' or 'bz2', or a `bool` or `NoneType`
        object representing the internal compression library for pickled `PackedPose` 
        objects and user-defined PyRosetta protocol `kwargs` objects. The default of
        `True` uses 'xz' for serialization if it's installed, otherwise uses 'zlib'
        for serialization.
        Default: True
    system_info: A `dict` or `NoneType` object specifying the system information
        required to reproduce the simulation. If `None` is provided, then PyRosettaCluster
        automatically detects the platform and returns this attribute as a dictionary
        {'sys.platform': `sys.platform`} (for example, {'sys.platform': 'linux'}).
        If a `dict` is provided, then validate that the 'sys.platform' key has a value
        equal to the current `sys.platform`, and log a warning message if not.
        Additional system information such as Amazon Machine Image (AMI) identifier
        and compute fleet instance type identifier may be stored in this dictionary,
        but is not validated. This information is stored in the simulation records for
        accounting.
        Default: None
    pyrosetta_build: A `str` or `NoneType` object specifying the PyRosetta build as
        output by `pyrosetta._build_signature()`. If `None` is provided, then PyRosettaCluster
        automatically detects the PyRosetta build and sets this attribute as the `str`.
        If a non-empty `str` is provided, then validate that the input PyRosetta build is
        equal to the active PyRosetta build, and raise an error if not. This ensures that
        reproduction simulations use an identical PyRosetta build from the original
        simulation. To bypass PyRosetta build validation with a warning message, an
        empty string ('') may be provided (but does not ensure reproducibility).
        Default: None
    sha1: A `str` or `NoneType` object specifying the git SHA1 hash string of the
        particular git commit being simulated. If a non-empty `str` object is provided,
        then it is validated to match the SHA1 hash string of the current HEAD,
        and then it is added to the simulation record for accounting. If an empty string
        is provided, then ensure that everything in the working directory is committed
        to the repository. If `None` is provided, then bypass SHA1 hash string
        validation and set this attribute to an empty string.
        Default: ""
    project_name: A `str` object specifying the project name of this simulation.
        This option just adds the user-provided project_name to the scorefile
        for accounting.
        Default: datetime.now().strftime("%Y.%m.%d.%H.%M.%S.%f") if not specified,
            else "PyRosettaCluster" if None
    simulation_name: A `str` object specifying the name of this simulation.
        This option just adds the user-provided `simulation_name` to the scorefile
        for accounting.
        Default: `project_name` if not specified, else "PyRosettaCluster" if None
    environment: A `NoneType` or `str` object specifying either the active conda/mamba environment
        YML file string, active uv project `requirements.txt` file string, or active pixi project
        `pixi.lock` file string. If a `NoneType` object is provided, then generate an environment file
        string for the active conda/mamba/uv/pixi environment and save it to the full simulation
        record. If a non-empty `str` object is provided, then validate it against the active
        conda/mamba/uv/pixi environment YML/requirements/lock file string and save it to the
        full simulation record. This ensures that reproduction simulations use an identical
        conda/mamba/uv/pixi environment to the original simulation. To bypass conda/mamba/uv/pixi
        environment validation with a warning message, an empty string ('') may be provided (but
        does not ensure reproducibility).
        Default: None
    output_path: A `str` object specifying the full path of the output directory
        (to be created if it doesn't exist) where the output results will be saved
        to disk.
        Default: "./outputs"
    output_init_file: A `str` object specifying the output ".init" file path that caches
        the 'input_packed_pose' keyword argument parameter upon PyRosettaCluster instantiation,
        and not including any output decoys, which is optionally used for exporting PyRosetta
        initialization files with output decoys by the `pyrosetta.distributed.cluster.export_init_file()`
        function after the simulation completes (see the 'output_decoy_types' keyword argument).
        If a `NoneType` object (or an empty `str` object ('')) is provided, or `dry_run=True`,
        then skip writing an output ".init" file upon PyRosettaCluster instantiation. If skipped,
        it is recommended to run `pyrosetta.dump_init_file()` before or after the simulation.
        If `compressed=True`, then the output file is further compressed by `bzip2`, and ".bz2"
        is appended to the filename.
        Default: `output_path`/`project_name`_`simulation_name`_pyrosetta.init
    output_decoy_types: An iterable of `str` objects representing the output decoy
        filetypes to save during the simulation. Available options are: ".pdb" for PDB
        files; ".pkl_pose" for pickled Pose files; ".b64_pose" for base64-encoded
        pickled Pose files; and ".init" for PyRosetta initialization files, each caching
        the host node PyRosetta initialization options (and input files, if any), the
        'input_packed_pose' keyword argument parameter (if any) and an output decoy.
        Because each ".init" file contains a copy of the PyRosetta initialization input files
        and input `PackedPose` object, unless these objects are relatively small in size
        or there are relatively few expected output decoys, then it is recommended to run
        `pyrosetta.distributed.cluster.export_init_file()` on only decoys of interest after the
        simulation completes without specifying ".init". If `compressed=True`, then each decoy
        output file is further compressed by `bzip2`, and ".bz2" is appended to the filename.
        Default: [".pdb",]
    output_scorefile_types: An iterable of `str` objects representing the output scorefile
        filetypes to save during the simulation. Available options are: ".json" for a
        JSON-encoded scorefile, and any filename extensions accepted by
        `pandas.DataFrame().to_pickle(compression="infer")` (including ".gz", ".bz2",
        and ".xz") for pickled `pandas.DataFrame` objects of scorefile data that can later
        be analyzed using `pyrosetta.distributed.cluster.io.secure_read_pickle(compression="infer")`.
        Note that in order to save pickled `pandas.DataFrame` objects, please ensure
        that `pyrosetta.secure_unpickle.add_secure_package("pandas")` has been first run.
        Default: [".json",]
    scorefile_name: A `str` object specifying the name of the output JSON-formatted
        scorefile, which must end in ".json". The scorefile location is always
        `output_path`/`scorefile_name`. If ".json" is not in the 'output_scorefile_types'
        keyword argument parameter, the JSON-formatted scorefile will not be output,
        but other scorefile types will get the same filename before the ".json" extension.
        Default: "scores.json"
    simulation_records_in_scorefile: A `bool` object specifying whether or not to
        write full simulation records to the scorefile. If `True`, then write
        full simulation records to the scorefile. This results in some redundant
        information on each line, allowing downstream reproduction of a decoy from
        the scorefile, but a larger scorefile. If `False`, then write
        curtailed simulation records to the scorefile. This results in minimally
        redundant information on each line, disallowing downstream reproduction
        of a decoy from the scorefile, but a smaller scorefile. If `False`, also
        write the active conda/mamba/uv/pixi environment to a file in the `output_path`
        keyword argument parameter. Full simulation records are always written to the
        output decoy files (the types of which are specified by the `output_decoy_types`
        keyword argument parameter), which can be used to reproduce any decoy without
        the scorefile.
        Default: False
    decoy_dir_name: A `str` object specifying the directory name where the
        output decoys will be saved. The directory location is always
        `output_path`/`decoy_dir_name`.
        Default: "decoys"
    logs_dir_name: A `str` object specifying the directory name where the
        output log files will be saved. The directory location is always
        `output_path`/`logs_dir_name`.
        Default: "logs"
    logging_level: A `str` object specifying the logging level of python tracer
        output to write to the log file of either "NOTSET", "DEBUG", "INFO",
        "WARNING", "ERROR", or "CRITICAL". The output log file is always written
        to `output_path`/`logs_dir_name`/`simulation_name`.log on disk.
        Default: "INFO"
    logging_address: A `str` object specifying the socket endpoint for sending and receiving
        log messages across a network, so log messages from user-provided PyRosetta
        protocols may be written to a single log file on the host node. The `str` object
        must take the format 'host:port' where 'host' is either an IP address, 'localhost',
        or Domain Name System (DNS)-accessible domain name, and the 'port' is a digit greater
        than or equal to 0. If the 'port' is '0', then the next free port is selected.
        Default: 'localhost:0' if `scheduler=None` or either the `client` or `clients`
            keyword argument parameters specify instances of `dask.distributed.LocalCluster`,
            otherwise '0.0.0.0:0'
    ignore_errors: A `bool` object specifying for PyRosettaCluster to ignore errors
        raised in the user-provided PyRosetta protocols. This comes in handy when
        well-defined errors are sparse and sporadic (such as rare Segmentation Faults),
        and the user would like PyRosettaCluster to run without raising the errors.
        Default: False
    timeout: A `float` or `int` object specifying how many seconds to wait between
        PyRosettaCluster checking-in on the running user-provided PyRosetta protocols.
        If each user-provided PyRosetta protocol is expected to run quickly, then
        0.1 seconds seems reasonable. If each user-provided PyRosetta protocol is
        expected to run slowly, then >1 second seems reasonable.
        Default: 0.5
    max_delay_time: A `float` or `int` object specifying the maximum number of seconds to 
        sleep before returning the result(s) from each user-provided PyRosetta protocol
        back to the client. If a dask worker returns the result(s) from a user-provided
        PyRosetta protocol too quickly, the dask scheduler needs to first register that
        the task is processing before it completes. In practice, in each user-provided
        PyRosetta protocol the runtime is subtracted from `max_delay_time`, and the dask
        worker sleeps for the remainder of the time, if any, before returning the result(s).
        It's recommended to set this option to at least 1 second, but longer times may
        be used as a safety throttle in cases of overwhelmed dask scheduler processes.
        Default: 3.0
    filter_results: A `bool` object specifying whether or not to filter out empty
        `PackedPose` objects between user-provided PyRosetta protocols. When a protocol
        returns or yields `NoneType`, PyRosettaCluster converts it to an empty `PackedPose`
        object that gets passed to the next protocol. If `True`, then filter out any empty
        `PackedPose` objects where there are no residues in the conformation as given by
        `Pose.empty()`, otherwise if `False` then continue to pass empty `PackedPose` objects
        to the next protocol. This is used for filtering out decoys mid-trajectory through
        user-provided PyRosetta protocols if protocols return or yield any `None`, empty
        `Pose`, or empty `PackedPose` objects.
        Default: True
    save_all: A `bool` object specifying whether or not to save all of the returned
        or yielded `Pose` and `PackedPose` objects from all user-provided
        PyRosetta protocols. This option may be used for checkpointing trajectories.
        To save arbitrary poses to disk, from within any user-provided PyRosetta
        protocol:
            `pose.dump_pdb(
                os.path.join(kwargs["PyRosettaCluster_output_path"], "checkpoint.pdb"))`
        Default: False
    dry_run: A `bool` object specifying whether or not to save '.pdb' files to
        disk. If `True`, then do not write '.pdb' or '.pdb.bz2' files to disk.
        Default: False
    security: A `bool` object or instance of `dask.distributed.Security()`, only having
        effect if `client=None` and `clients=None`, that is passed to 'dask' if using
        `scheduler=None` or passed to 'dask-jobqueue' if using `scheduler="slurm"` or
        `scheduler="sge"`. If `True` is provided, then invoke the 'cryptography' package
        to generate a `Security.temporary()` object through 'dask' or 'dask-jobqueue'. See
        https://distributed.dask.org/en/latest/tls.html#distributed.security.Security.temporary
        for more information. If a dask `Security()` object is provided, then pass it to
        dask with `scheduler=None`, or pass it to 'dask-jobqueue' (where 'shared_temp_directory'
        is set to the `output_path` keyword argument parameter) with `scheduler="slurm"` or
        `scheduler="sge"`. If `False` is provided, then security is disabled regardless
        of the `scheduler` keyword argument parameter (which is not recommended for remote
        clusters unless using a firewall). If `None` is provided, then `True` is used by
        default. In order to generate a `dask.distributed.Security()` object with OpenSSL,
        the `pyrosetta.distributed.cluster.generate_dask_tls_security()` function may also
        be used (see docstring for more information) instead of the 'cryptography' package.
        Default: `False` if `scheduler=None`, otherwise `True`
    max_nonce: An `int` object greater than or equal to 1 specifying the maximum number of
        nonces to cache per process if dask security is disabled while using remote clusters,
        which protects against replay attacks. If nonce caching is in use, each process
        (including the host node process and all dask worker processes) cache nonces upon
        communication exchange over the network, which can increase memory usage in each
        process. A rough estimate of additional memory usage is ~0.2 KB per task
        per user-provided PyRosetta protocol per process. For example, submitting
        1000 tasks with 2 user-provided PyRosetta protocols adds ~0.2 KB/task/protocol
        * 1000 tasks * 2 protocols = ~0.4 MB of memory per processs. If memory usage
        per process permits, it is recommended to set this parameter to at least the
        number of tasks times the number of protocols submitted, so that every nonce
        from every communication exchange over the network gets cached.
        Default: 4096
    cooldown_time: A `float` or `int` object specifying how many seconds to sleep after the
        simulation is complete to allow loggers to flush. For very slow network filesystems,
        2.0 or more seconds may be reasonable.
        Default: 0.5
    norm_task_options: A `bool` object specifying whether or not to normalize the task
        'options' and 'extra_options' values after PyRosetta initialization on the remote
        compute cluster. If `True`, then this enables more facile simulation reproduction
        by the use of the `ProtocolSettingsMetric` SimpleMetric to normalize the PyRosetta
        initialization options and by relativization of any input files and directory paths
        to the current working directory from which the task is running.
        Default: True
    author: An optional `str` object specifying the author(s) of the simulation that is
        written to the full simulation records and the PyRosetta initialization '.init' file.
        Default: ""
    email: An optional `str` object specifying the email address(es) of the author(s) of
        the simulation that is written to the full simulation records and the PyRosetta
        initialization '.init' file.
        Default: ""
    license: An optional `str` object specifying the license of the output data of the
        simulation that is written to the full simulation records and the PyRosetta
        initialization '.init' file (e.g., "ODC-ODbL", "CC BY-ND", "CDLA Permissive-2.0", etc.).
        Default: ""

Returns:
    A PyRosettaCluster instance.
"""
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import attr
    import distributed
    import toolz
    from dask.distributed import Client, Future, Security, as_completed
    from distributed.scheduler import KilledWorker
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.core' requires the "
        + "third-party packages 'attrs', 'dask', 'distributed', and 'toolz' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/attrs/\n"
        + "https://pypi.org/project/dask/\n"
        + "https://pypi.org/project/distributed/\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import logging
import os
import uuid

from datetime import datetime
from pyrosetta.distributed.cluster.base import TaskBase, _get_residue_type_set
from pyrosetta.distributed.cluster.config import get_environment_manager
from pyrosetta.distributed.cluster.converters import (
    is_empty as _is_empty,
    _maybe_issue_environment_warnings,
    _parse_decoy_ids,
    _parse_filter_results,
    _parse_environment,
    _parse_input_packed_pose,
    _parse_logging_address,
    _parse_norm_task_options,
    _parse_output_decoy_types,
    _parse_output_scorefile_types,
    _parse_pyrosetta_build,
    _parse_scratch_dir,
    _parse_seeds,
    _parse_sha1,
    _parse_system_info,
    _parse_tasks,
    _parse_yield_results,
)
from pyrosetta.distributed.cluster.hkdf import derive_task_key
from pyrosetta.distributed.cluster.initialization import _get_pyrosetta_init_args, _maybe_init_client
from pyrosetta.distributed.cluster.io import IO
from pyrosetta.distributed.cluster.logging_support import LoggingSupport, MaskedBytes
from pyrosetta.distributed.cluster.multiprocessing import user_spawn_thread
from pyrosetta.distributed.cluster.security import SecurityIO
from pyrosetta.distributed.cluster.serialization import NonceCache, Serialization
from pyrosetta.distributed.cluster.utilities import SchedulerManager
from pyrosetta.distributed.cluster.validators import (
    _validate_dir,
    _validate_dirs,
    _validate_float,
    _validate_int,
    _validate_logging_address,
    _validate_min_len,
    _validate_output_init_file,
    _validate_scorefile_name,
)
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import (
    Any,
    Dict,
    Generator,
    List,
    NoReturn,
    Optional,
    Tuple,
    TypeVar,
    Union,
)


G = TypeVar("G")


@attr.s(kw_only=True, slots=True, frozen=False)
class PyRosettaCluster(IO[G], LoggingSupport[G], SchedulerManager[G], SecurityIO[G], TaskBase[G]):
    tasks = attr.ib(
        type=list,
        default=[{}],
        validator=attr.validators.deep_iterable(
            member_validator=attr.validators.instance_of(dict),
            iterable_validator=attr.validators.instance_of(list),
        ),
        converter=_parse_tasks,
    )
    nstruct = attr.ib(
        type=int,
        default=1,
        validator=[_validate_int, attr.validators.instance_of(int)],
        converter=attr.converters.default_if_none(default=1),
    )
    tasks_size = attr.ib(
        type=int,
        default=attr.Factory(
            lambda self: toolz.itertoolz.count(self.tasks) * self.nstruct,
            takes_self=True,
        ),
        init=False,
        validator=attr.validators.instance_of(int),
    )
    input_packed_pose = attr.ib(
        type=PackedPose,
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(PackedPose)),
        converter=_parse_input_packed_pose,
    )
    seeds = attr.ib(
        type=list,
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(list)),
        converter=attr.converters.optional(_parse_seeds),
    )
    decoy_ids = attr.ib(
        type=list,
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(list)),
        converter=attr.converters.optional(_parse_decoy_ids),
    )
    client = attr.ib(
        type=Optional[distributed.client.Client],
        default=None,
        validator=attr.validators.optional(
            attr.validators.instance_of(distributed.client.Client)
        ),
    )
    clients = attr.ib(
        type=Optional[List[distributed.client.Client]],
        default=None,
        validator=[
            attr.validators.optional(
                attr.validators.deep_iterable(
                    member_validator=attr.validators.instance_of(distributed.client.Client),
                    iterable_validator=attr.validators.instance_of((tuple, list)),
                ),
            ),
            _validate_min_len,
        ],
    )
    scheduler = attr.ib(
        type=str,
        default=None,
        validator=attr.validators.optional(
            [attr.validators.in_(["sge", "slurm"]), attr.validators.instance_of(str),],
        ),
    )
    cores = attr.ib(
        type=int,
        default=1,
        validator=[_validate_int, attr.validators.instance_of(int)],
        converter=attr.converters.default_if_none(default=1),
    )
    processes = attr.ib(
        type=int,
        default=1,
        validator=[_validate_int, attr.validators.instance_of(int)],
        converter=attr.converters.default_if_none(default=1),
    )
    memory = attr.ib(
        type=str,
        default="4g",
        validator=attr.validators.optional(attr.validators.instance_of(str)),
        converter=attr.converters.default_if_none(default="4g"),
    )
    scratch_dir = attr.ib(
        type=str,
        default=None,
        validator=[_validate_dir, attr.validators.instance_of(str)],
        converter=_parse_scratch_dir,
    )
    adapt_threshold = attr.ib(
        type=int,
        default=1000,
        init=False,
        validator=[_validate_int, attr.validators.instance_of(int)],
    )
    min_workers = attr.ib(
        type=int,
        default=1,
        validator=[_validate_int, attr.validators.instance_of(int)],
        converter=attr.converters.default_if_none(default=1),
    )
    max_workers = attr.ib(
        type=int,
        default=attr.Factory(
            lambda self: self.adapt_threshold if (self.tasks_size < self.adapt_threshold) else self.tasks_size,
            takes_self=True,
        ),
        validator=[_validate_int, attr.validators.instance_of(int)],
        converter=attr.converters.default_if_none(default=1000),
    )
    dashboard_address = attr.ib(
        type=str,
        default=":8787",
        validator=attr.validators.instance_of(str),
        converter=attr.converters.default_if_none(default=":8787"),
    )
    project_name = attr.ib(
        type=str,
        default=datetime.now().strftime("%Y.%m.%d.%H.%M.%S.%f"),
        validator=attr.validators.optional(attr.validators.instance_of(str)),
        converter=attr.converters.default_if_none(default="PyRosettaCluster"),
    )
    simulation_name = attr.ib(
        type=str,
        default=attr.Factory(lambda self: self.project_name, takes_self=True),
        validator=attr.validators.instance_of(str),
        converter=attr.converters.default_if_none(default="PyRosettaCluster"),
    )
    output_path = attr.ib(
        type=str,
        default=os.path.abspath(os.path.join(os.getcwd(), "outputs")),
        validator=[_validate_dirs, attr.validators.instance_of(str)],
        converter=attr.converters.default_if_none(
            default=os.path.abspath(os.path.join(os.getcwd(), "outputs"))
        ),
    )
    output_decoy_types = attr.ib(
        type=List[str],
        default=None,
        validator=[
            attr.validators.deep_iterable(
                member_validator=attr.validators.instance_of(str),
                iterable_validator=attr.validators.instance_of(list),
            ),
            _validate_min_len,
        ],
        converter=_parse_output_decoy_types,
    )
    output_scorefile_types = attr.ib(
        type=List[str],
        default=None,
        validator=[
            attr.validators.deep_iterable(
                member_validator=attr.validators.instance_of(str),
                iterable_validator=attr.validators.instance_of(list),
            ),
            _validate_min_len,
        ],
        converter=_parse_output_scorefile_types,
    )
    scorefile_name = attr.ib(
        type=str,
        default="scores.json",
        validator=[attr.validators.instance_of(str), _validate_scorefile_name],
        converter=attr.converters.default_if_none(default="scores.json"),
    )
    scorefile_path = attr.ib(
        type=str,
        default=attr.Factory(
            lambda self: os.path.abspath(
                os.path.join(self.output_path, self.scorefile_name)
            ),
            takes_self=True,
        ),
        init=False,
        validator=attr.validators.instance_of(str),
    )
    simulation_records_in_scorefile = attr.ib(
        type=bool,
        default=False,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    decoy_dir_name = attr.ib(
        type=str,
        default="decoys",
        validator=attr.validators.instance_of(str),
        converter=attr.converters.default_if_none(default="decoys"),
    )
    decoy_path = attr.ib(
        type=str,
        default=attr.Factory(
            lambda self: os.path.abspath(
                os.path.join(self.output_path, self.decoy_dir_name)
            ),
            takes_self=True,
        ),
        init=False,
        validator=attr.validators.instance_of(str),
    )
    logs_dir_name = attr.ib(
        type=str,
        default="logs",
        validator=attr.validators.instance_of(str),
        converter=attr.converters.default_if_none(default="logs"),
    )
    logs_path = attr.ib(
        type=str,
        default=attr.Factory(
            lambda self: os.path.abspath(
                os.path.join(self.output_path, self.logs_dir_name)
            ),
            takes_self=True,
        ),
        init=False,
        validator=attr.validators.instance_of(str),
    )
    logging_level = attr.ib(
        type=str,
        default="INFO",
        validator=[
            attr.validators.in_(
                ["NOTSET", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
            ),
            attr.validators.instance_of(str),
        ],
        converter=attr.converters.default_if_none(default="NOTSET"),
    )
    logging_file = attr.ib(
        type=str,
        default=attr.Factory(
            lambda self: os.path.join(
                self.logs_path,
                "_".join(
                    [
                        self.project_name.replace(" ", "-"),
                        self.simulation_name.replace(" ", "-") + ".log",
                    ]
                ),
            ),
            takes_self=True,
        ),
        init=False,
        validator=attr.validators.instance_of(str),
    )
    logging_address = attr.ib(
        type=str,
        default=attr.Factory(_parse_logging_address, takes_self=True),
        validator=[_validate_logging_address, attr.validators.instance_of(str)],
    )
    compressed = attr.ib(
        type=bool,
        default=True,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=True),
    )
    compression = attr.ib(
        type=Optional[Union[str, bool]],
        default=True,
        validator=attr.validators.optional(attr.validators.instance_of((bool, str))),
    )
    sha1 = attr.ib(
        type=str,
        default="",
        validator=attr.validators.instance_of(str),
        converter=_parse_sha1,
    )
    ignore_errors = attr.ib(
        type=bool,
        default=False,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    timeout = attr.ib(
        type=float,
        default=0.5,
        validator=[_validate_float, attr.validators.instance_of((float, int))],
        converter=attr.converters.default_if_none(default=0.5),
    )
    max_delay_time = attr.ib(
        type=float,
        default=3.0,
        validator=[_validate_float, attr.validators.instance_of((float, int))],
        converter=attr.converters.default_if_none(default=3.0),
    )
    filter_results = attr.ib(
        type=bool,
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=_parse_filter_results,
    )
    save_all = attr.ib(
        type=bool,
        default=False,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(default=False),
    )
    dry_run = attr.ib(
        type=bool,
        default=False,
        validator=attr.validators.instance_of(bool),
        converter=attr.converters.default_if_none(False),
    )
    norm_task_options = attr.ib(
        type=bool,
        default=None,
        validator=attr.validators.instance_of(bool),
        converter=_parse_norm_task_options,
    )
    yield_results = attr.ib(
        type=bool,
        init=False,
        default=False,
        validator=attr.validators.instance_of(bool),
        converter=_parse_yield_results,
    )
    cooldown_time = attr.ib(
        type=float,
        default=0.5,
        validator=[_validate_float, attr.validators.instance_of((float, int))],
        converter=attr.converters.default_if_none(default=0.5),
    )
    protocols_key = attr.ib(
        type=str,
        default="PyRosettaCluster_protocols_container",
        init=False,
        validator=attr.validators.instance_of(str),
    )
    system_info = attr.ib(
        type=dict,
        default=None,
        validator=attr.validators.instance_of(dict),
        converter=_parse_system_info,
    )
    pyrosetta_build = attr.ib(
        type=str,
        default=None,
        validator=attr.validators.instance_of(str),
        converter=_parse_pyrosetta_build,
    )
    security = attr.ib(
        type=Union[bool, Security],
        default=attr.Factory(
            lambda self: bool(self.scheduler),
            takes_self=True,
        ),
        validator=attr.validators.instance_of((bool, Security)),
        converter=attr.converters.default_if_none(default=True),
    )
    instance_id = attr.ib(
        type=str,
        default=attr.Factory(
            lambda: f"PyRosettaCluster_instance_id_{uuid.uuid4().hex}",
            takes_self=False,
        ),
        init=False,
        validator=attr.validators.instance_of(str),
    )
    max_nonce = attr.ib(
        type=int,
        default=4096,
        validator=[attr.validators.instance_of(int), _validate_int],
    )
    environment = attr.ib(
        type=str,
        default=None,
        validator=attr.validators.instance_of(str),
        converter=_parse_environment,
    )
    author = attr.ib(
        type=str,
        default=None,
        validator=attr.validators.instance_of(str),
        converter=attr.converters.default_if_none(""),
    )
    email = attr.ib(
        type=str,
        default=None,
        validator=attr.validators.instance_of(str),
        converter=attr.converters.default_if_none(""),
    )
    license = attr.ib(
        type=str,
        default=None,
        validator=attr.validators.instance_of(str),
        converter=attr.converters.default_if_none(""),
    )
    output_init_file = attr.ib(
        type=str,
        default=attr.Factory(
            lambda self: os.path.join(
                self.output_path,
                "_".join(
                    [
                        self.project_name.replace(" ", "-"),
                        self.simulation_name.replace(" ", "-"),
                        "pyrosetta.init",
                    ]
                ),
            ),
            takes_self=True,
        ),
        validator=[attr.validators.instance_of(str), _validate_output_init_file],
        converter=attr.converters.default_if_none(""),
    )
    environment_manager = attr.ib(
        type=str,
        default=attr.Factory(
            get_environment_manager,
            takes_self=False,
        ),
        init=False,
        validator=attr.validators.instance_of(str),
    )
    environment_file = attr.ib(
        type=str,
        default=attr.Factory(
            lambda self: os.path.join(
                self.output_path,
                "_".join(
                    [
                        self.project_name.replace(" ", "-"),
                        self.simulation_name.replace(" ", "-"),
                        (
                            "pixi.lock"
                            if self.environment_manager == "pixi"
                            else "requirements.txt"
                            if self.environment_manager == "uv"
                            else "environment.yml"
                        )
                    ]
                ),
            ),
            takes_self=True,
        ),
        init=False,
        validator=attr.validators.instance_of(str),
    )
    pyrosetta_init_args = attr.ib(
        type=list,
        default=attr.Factory(_get_pyrosetta_init_args),
        init=False,
        validator=attr.validators.deep_iterable(
            member_validator=attr.validators.instance_of(str),
            iterable_validator=attr.validators.instance_of(list),
        ),
    )

    def __attrs_pre_init__(self) -> None:
        _maybe_issue_environment_warnings()

    def __attrs_post_init__(self) -> None:
        _maybe_init_client()
        self._setup_logger()
        self._cache_toml()
        self._write_environment_file(self.environment_file)
        self._write_init_file()
        self.serializer = Serialization(compression=self.compression)
        self.clients_dict = self._setup_clients_dict()
        self.with_nonce = self._setup_with_nonce()
        self.serializer = Serialization(
            instance_id=self.instance_id,
            compression=self.compression,
            with_nonce=self.with_nonce,
        )
        self.nonce_cache = NonceCache(
            instance_id=self.instance_id,
            max_nonce=self.max_nonce,
        ) if self.with_nonce else None

    def _create_future(
        self,
        client: Client,
        protocol_name: str,
        compressed_protocol: bytes,
        compressed_packed_pose: bytes,
        compressed_kwargs: bytes,
        pyrosetta_init_kwargs: Dict[str, Any],
        extra_args: Dict[str, Any],
        passkey: bytes,
        resource: Optional[Dict[Any, Any]],
        priority: Optional[int],
        retry: Optional[int],
    ) -> Future:
        """Scatter data and return submitted 'user_spawn_thread' future."""
        task_id = uuid.uuid4().hex
        masked_key = MaskedBytes(derive_task_key(passkey, task_id))
        scatter = client.scatter(
            (
                protocol_name,
                compressed_protocol,
                compressed_packed_pose,
                compressed_kwargs,
                pyrosetta_init_kwargs,
                repr(client),
                extra_args,
                masked_key,
                task_id,
            ),
            broadcast=False,
            hash=False,
        )
        submit_kwargs = {"pure": False}
        # Omit resources keyword argument for distributed versions <2.1.0
        # or use default if user specifies `resources=None` in distributed versions >=2.1.0
        if resource is not None:
            submit_kwargs["resources"] = resource
        # Omit priority keyword argument for distributed versions <1.21.0
        # or use default if user specifies `priorities=None` in distributed versions >=1.21.0
        if priority is not None:
            submit_kwargs["priority"] = priority
        # Omit retries keyword argument for distributed versions <1.20.0
        # or use default if user specifies `retries=None` for distributed versions >=1.20.0
        if retry is not None:
            submit_kwargs["retries"] = retry

        return client.submit(user_spawn_thread, *scatter, **submit_kwargs)

    def _run(
        self,
        *args: Any,
        protocols: Any = None,
        clients_indices: Any = None,
        resources: Any = None,
        priorities: Any = None,
        retries: Any = None,
    ) -> Union[NoReturn, Generator[Tuple[PackedPose, Dict[Any, Any]], None, None]]:
        """
        Run user-provided PyRosetta protocols on a local or remote compute cluster using
        the user-customized PyRosettaCluster instance. Either arguments or the 'protocols'
        keyword argument is required. If both are provided, then the 'protocols' keyword
        argument gets concatenated after the input arguments.

        Examples:
            PyRosettaCluster().distribute(protocol_1)
            PyRosettaCluster().distribute(protocols=protocol_1)
            PyRosettaCluster().distribute(protocol_1, protocol_2, protocol_3)
            PyRosettaCluster().distribute(protocols=(protocol_1, protocol_2, protocol_3))
            PyRosettaCluster().distribute(protocol_1, protocol_2, protocols=[protocol_3, protocol_4])

            # Run `protocol_1` on `client_1`,
            # then `protocol_2` on `client_2`, 
            # then `protocol_3` on `client_1`,
            # then `protocol_4` on `client_2`:
            PyRosettaCluster(clients=[client_1, client_2]).distribute(
                protocols=[protocol_1, protocol_2, protocol_3, protocol_4],
                clients_indices=[0, 1, 0, 1],
            )

            # Run `protocol_1` on `client_2`,
            # then `protocol_2` on `client_3`,
            # then `protocol_3` on `client_1`:
            PyRosettaCluster(clients=[client_1, client_2, client_3]).distribute(
                protocols=[protocol_1, protocol_2, protocol_3],
                clients_indices=[1, 2, 0],
            )

            # Run `protocol_1` on `client_1` with dask worker resource constraints "GPU=2",
            # then `protocol_2` on `client_1` with dask worker resource constraints "MEMORY=100e9",
            # then `protocol_3` on `client_1` without dask worker resource constraints:
            PyRosettaCluster(client=client_1).distribute(
                protocols=[protocol_1, protocol_2, protocol_3],
                resources=[{"GPU": 2}, {"MEMORY": 100e9}, None],
            )

            # Run `protocol_1` on `client_1` with dask worker resource constraints "GPU=2",
            # then `protocol_2` on `client_2` with dask worker resource constraints "MEMORY=100e9":
            PyRosettaCluster(clients=[client_1, client_2]).distribute(
                protocols=[protocol_1, protocol_2],
                clients_indices=[0, 1],
                resources=[{"GPU": 2}, {"MEMORY": 100e9}],
            )

            # Run protocols with depth-first task execution:
            PyRosettaCluster().distribute(
                protocols=[protocol_1, protocol_2, protocol_3, protocol_4],
                priorities=[0, 10, 20, 30],
            )

            # Run protocols with up to three retries per failed task during `protocol_3` and `protocol_4`:
            PyRosettaCluster().distribute(
                protocols=[protocol_1, protocol_2, protocol_3, protocol_4],
                retries=[0, 0, 3, 3],
            )

        Args:
            *args: Optional instances of type `types.GeneratorType` or `types.FunctionType`,
                in the order of protocols to be executed.
            protocols: An optional iterable of extra callable PyRosetta protocols,
                i.e. an iterable of objects of `types.GeneratorType` and/or
                `types.FunctionType` types; or a single instance of type
                `types.GeneratorType` or `types.FunctionType`.
                Default: None
            clients_indices: An optional `list` or `tuple` object of `int` objects, where each `int` object represents
                a zero-based index corresponding to the initialized dask `distributed.client.Client` object(s) passed 
                to the `PyRosettaCluster(clients=...)` keyword argument parameter. If not `None`, then the length of the
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
                applied, the protocols will not run. See https://distributed.dask.org/en/stable/resources.html for more
                information.
                Default: None
            priorities: An optional `list` or `tuple` of `int` objects, where each `int` object sets the dask scheduler
                priority for the corresponding user-defined PyRosetta protocol (i.e., indexed the same as `client_indices`).
                If `None`, then no explicit priorities are set. If not `None`, then the length of the `priorities` object
                must equal the number of protocols passed to the `PyRosettaCluster().distribute` method, and each `int`
                value determines the dask scheduler priority for that protocol's received tasks.
                Breadth-first task execution (default):
                    When all user-defined PyRosetta protocols have an identical priority (e.g., `[0] * len(protocols)` or
                    `None`), then all tasks enter the dask scheduler's queue with equal priority. Under equal priority, dask
                    mainly schedules tasks in a first-in, first-out manner. When dask worker resources are saturated, this
                    causes all tasks submitted to upstream protocols to run to completion before tasks are scheduled to run
                    downstream protocols, producing a breadth-first task execution behavior across user-defined PyRosetta
                    protocols.
                Depth-first task execution:
                    To allow tasks to run through all user-defined PyRosetta protocols before all tasks submitted to
                    upstream protocols complete, assign increasing priorities to downstream protocols (e.g.,
                    `list(range(0, len(protocols) * 10, 10))`). Once a task completes an upstream protocol, it is submitted
                    to the next downstream protocol with a higher priority than tasks still queued for upstream protocols,
                    so tasks may run through all user-defined PyRosetta protocols to completion as dask worker resources
                    become available. This produces a depth-first task execution behavior across user-defined PyRosetta
                    protocols when dask worker resources are saturated.
                See https://distributed.dask.org/en/stable/priority.html for more information.
                Default: None
            retries: An optional `list` or `tuple` of `int` objects, where each `int` object (≥0) sets the number of allowed
                automatic retries of each failed task that was submitted to the corresponding user-provided PyRosetta protocol
                (i.e., indexed the same as `client_indices`). If an `int` object (≥0) is provided, then apply that number of
                allowed automatic retries to all user-provided PyRosetta protocols. If `None`, then no explicit retries are
                allowed. If not `None` and not an `int` object, then the length of the `retries` parameter must equal the number
                of protocols passed to the `PyRosettaCluster().distribute` method, and each `int` value determines the number
                of automatic retries the dask scheduler allows for that protocol's failed tasks. Allowing retries of failed tasks
                may be useful if remote compute resources are subject to preemption (e.g., cloud spot instances or backfill
                queues). Note that retries are only appropriate for user-provided PyRosetta protocols that are side effect-free
                upon preemption, in which tasks can be restarted without producing inconsistent external states if preempted midway
                through the protocol. Also note that if `PyRosettaCluster(ignore_errors=True)` is used, then protocols failing due
                to standard Python exceptions or Rosetta segmentation faults will still be considered successes, and this
                keyword argument parameter has no effect on them since these protocol errors are ignored. However, if a compute
                resource executing tasks is reclaimed midway through a protocol, then the dask scheduler registers those tasks
                as incomplete, and they may be retried a certain number of times based on this keyword argument parameter.
                See https://distributed.dask.org/en/latest/scheduling-state.html#task-state for more information.
                Default: None
        """
        yield_results = _parse_yield_results(self.yield_results)
        clients, cluster, adaptive = self._setup_clients_cluster_adaptive()
        self._setup_task_security_plugin(clients)
        socket_listener_address, passkey = self._setup_socket_listener(clients)
        compressed_input_packed_pose = self.serializer.compress_packed_pose(self.input_packed_pose)
        resources = self._parse_resources(resources)
        priorities = self._parse_priorities(priorities)
        retries = self._parse_retries(retries)
        protocols, protocol, seed, clients_index, resource, priority, retry = self._setup_protocols_protocol_seed(
            args, protocols, clients_indices, resources, priorities, retries
        )
        protocol_name = protocol.__name__
        compressed_protocol = self.serializer.compress_object(protocol)
        client_residue_type_set = _get_residue_type_set()
        extra_args = dict(
            decoy_ids=self.decoy_ids,
            protocols_key=self.protocols_key,
            timeout=self.timeout,
            ignore_errors=self.ignore_errors,
            datetime_format=self.DATETIME_FORMAT,
            instance_id=self.instance_id,
            compression=self.compression,
            with_nonce=self.with_nonce,
            norm_task_options=self.norm_task_options,
            max_delay_time=self.max_delay_time,
            logging_level=self.logging_level,
            socket_listener_address=socket_listener_address,
            client_residue_type_set=client_residue_type_set,
        )
        seq = as_completed(
            [
                self._create_future(
                    clients[clients_index],
                    protocol_name,
                    compressed_protocol,
                    compressed_input_packed_pose,
                    compressed_kwargs,
                    pyrosetta_init_kwargs,
                    extra_args,
                    passkey,
                    resource,
                    priority,
                    retry,
                )
                for _ in range(self.nstruct)
                for compressed_kwargs, pyrosetta_init_kwargs in (
                    self._setup_initial_kwargs(protocols, seed, task_kwargs) for task_kwargs in self.tasks
                )
            ]
        )
        for i, future in enumerate(seq, start=1):
            try:
                results = future.result()
            except KilledWorker as ex:
                logging.error(ex)
                continue
            logging.info(
                "Percent Complete = {0:0.5f} %".format((i / self.tasks_size) * 100.0)
            )
            for compressed_packed_pose, compressed_kwargs in results:
                if self.with_nonce:
                    self.nonce_cache._cache_nonce(compressed_kwargs)
                kwargs = self.serializer.decompress_kwargs(compressed_kwargs)
                if not kwargs[self.protocols_key]:
                    if yield_results:
                        yield self.serializer.decompress_packed_pose(compressed_packed_pose), self.serializer.deepcopy_kwargs(kwargs)
                    self._save_results(compressed_packed_pose, kwargs)
                else:
                    if self.save_all:
                        if yield_results:
                            yield self.serializer.decompress_packed_pose(compressed_packed_pose), self.serializer.deepcopy_kwargs(kwargs)
                        self._save_results(
                            compressed_packed_pose,
                            self.serializer.deepcopy_kwargs(kwargs),
                        )
                    if self.filter_results and _is_empty(self.serializer.decompress_packed_pose(compressed_packed_pose)):
                        continue
                    compressed_kwargs, pyrosetta_init_kwargs, protocol, clients_index, resource, priority, retry = self._setup_kwargs(
                        kwargs, clients_indices, resources, priorities, retries
                    )
                    protocol_name = protocol.__name__
                    compressed_protocol = self.serializer.compress_object(protocol)
                    seq.add(
                        self._create_future(
                            clients[clients_index],
                            protocol_name,
                            compressed_protocol,
                            compressed_packed_pose,
                            compressed_kwargs,
                            pyrosetta_init_kwargs,
                            extra_args,
                            passkey,
                            resource,
                            priority,
                            retry,
                        )
                    )
                    self.tasks_size += 1
                    self._maybe_adapt(adaptive)

        self._close_socket_listener(clients)
        self._maybe_teardown(clients, cluster)
        self._close_logger()

    def generate(
        self,
        *args: Any,
        protocols: Any = None,
        clients_indices: Any = None,
        resources: Any = None,
        priorities: Any = None,
        retries: Any = None,
    ) -> Union[NoReturn, Generator[Tuple[PackedPose, Dict[Any, Any]], None, None]]:
        if self.sha1 != "":
            logging.warning(
                "Use of the `PyRosettaCluster.generate` method for reproducible simulations is not supported! "
                + "PyRosettaCluster reproduces decoys from the output files written to disk. Subsequent code run "
                + "on these results is not being saved by PyRosettaCluster. Use the `PyRosettaCluster.distribute` "
                + "method for reproducible simulations. To silence this warning and continue without using version "
                + "control, set the `sha1` instance attribute of PyRosettaCluster to `None`."
            )
        self.yield_results = True
        for result in self._run(
            *args,
            protocols=protocols,
            clients_indices=clients_indices,
            resources=resources,
            priorities=priorities,
            retries=retries,
        ):
            yield result

    def distribute(
        self,
        *args: Any,
        protocols: Any = None,
        clients_indices: Any = None,
        resources: Any = None,
        priorities: Any = None,
        retries: Any = None,
    ) -> Optional[NoReturn]:
        self.yield_results = False
        for _ in self._run(
            *args,
            protocols=protocols,
            clients_indices=clients_indices,
            resources=resources,
            priorities=priorities,
            retries=retries,
        ):
            pass

    distribute.__doc__ = _run.__doc__ + """
        Returns:
            None
        """

    generate.__doc__ = _run.__doc__ + """
        Extra information:

        The `PyRosettaCluster.generate` method may be used for developing PyRosetta protocols on a local or
        remote compute cluster and optionally post-processing or visualizing output `PackedPose` objects in memory.
        Importantly, subsequent code run on the yielded results is not captured by PyRosettaCluster, and so use
        of this method does not ensure reproducibility of the simulation. Use the `PyRosettaCluster.distribute`
        method for reproducible simulations.

        Each yielded result is a `tuple` object with a `PackedPose` object as the first element and a `dict`
        object as the second element. The `PackedPose` object represents a returned or yielded `PackedPose`
        (or `Pose` or `NoneType`) object from the most recently run user-provided PyRosetta protocol. The `dict`
        object represents the optionally returned or yielded user-defined PyRosetta protocol `kwargs` dictionary
        object from the same most recently run user-provided PyRosetta protocol (see 'protocols' argument). If
        `PyRosettaCluster(save_all=True)`, results are yielded after each user-provided PyRosetta protocol,
        otherwise results are yielded after the final user-defined PyRosetta protocol. Results are yielded in the
        order in which they arrive back to the client(s) from the distributed cluster (which may differ from the
        order that tasks are submitted, due to tasks running asynchronously). If `PyRosettaCluster(dry_run=True)`,
        results are still yielded but '.pdb' or '.pdb.bz2' files are not saved to disk.
        See https://docs.dask.org/en/latest/futures.html#distributed.as_completed for more information.

        Extra examples:

            # Iterate over results in real-time as they are yielded from the cluster:
            for packed_pose, kwargs in PyRosettaCluster().generate(protocols):
                ...

            # Iterate over submissions to the same client:
            client = Client()
            for packed_pose, kwargs in PyRosettaCluster(client=client).generate(protocols):
                # Post-process results on host node asynchronously from results generation
                prc = PyRosettaCluster(
                    input_packed_pose=packed_pose,
                    client=client,
                    logs_dir_name=f"logs_{uuid.uuid4().hex}", # Make sure to write new log files
                )
                for packed_pose, kwargs in prc.generate(other_protocols):
                    ...

            # Iterate over multiple clients:
            client_1 = Client()
            client_2 = Client()
            for packed_pose, kwargs in PyRosettaCluster(client=client_1).generate(protocols):
                # Post-process results on host node asynchronously from results generation
                prc = PyRosettaCluster(
                    input_packed_pose=packed_pose,
                    client=client_2,
                    logs_dir_name=f"logs_{uuid.uuid4().hex}", # Make sure to write new log files
                )
                for packed_pose, kwargs in prc.generate(other_protocols):
                    ...

            # Using multiple `dask.distributed.as_completed` iterators on the host node creates additional overhead.
            # If post-processing on the host node is not required between user-provided PyRosetta protocols,
            # the preferred method is to distribute protocols within a single `PyRosettaCluster().generate()`
            # method call using the `clients_indices` keyword argument:
            prc_generate = PyRosettaCluster(clients=[client_1, client_2]).generate(
                protocols=[protocol_1, protocol_2],
                clients_indices=[0, 1],
                resources=[{"GPU": 1}, {"CPU": 1}],
            )
            for packed_pose, kwargs in prc_generate:
                # Post-process results on host node asynchronously from results generation

        Yields:
            (PackedPose, dict) tuples from the most recently run user-provided PyRosetta protocol if
            `PyRosettaCluster(save_all=True)` otherwise from the final user-defined PyRosetta protocol.
        """

PyRosettaCluster.__doc__ = __doc__
