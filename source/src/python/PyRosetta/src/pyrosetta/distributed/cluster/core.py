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
        '.pdb' files with bzip2, resulting in '.pdb.bz2' files.
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
        output by `pyrosetta._version_string()`. If `None` is provided, then PyRosettaCluster
        automatically detects the PyRosetta build and sets this attribute as the `str`.
        If a `str` is provided, then validate that the input PyRosetta build is equal
        to the active PyRosetta build, and log a warning message if not.
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
    environment: A `NoneType` or `str` object specifying the active conda environment
        YML file string. If a `NoneType` object is provided, then generate a YML file
        string for the active conda environment and save it to the full simulation
        record. If a `str` object is provided, then validate it against the active
        conda environment YML file string and save it to the full simulation record.
        Default: None
    output_path: A `str` object specifying the full path of the output directory
        (to be created if it doesn't exist) where the output results will be saved
        to disk.
        Default: "./outputs"
    scorefile_name: A `str` object specifying the name of the output JSON-formatted
        scorefile. The scorefile location is always `output_path`/`scorefile_name`.
        Default: "scores.json"
    simulation_records_in_scorefile: A `bool` object specifying whether or not to
        write full simulation records to the scorefile. If `True`, then write
        full simulation records to the scorefile. This results in some redundant
        information on each line, allowing downstream reproduction of a decoy from
        the scorefile, but a larger scorefile. If `False`, then write
        curtailed simulation records to the scorefile. This results in minimally
        redundant information on each line, disallowing downstream reproduction
        of a decoy from the scorefile, but a smaller scorefile. If `False`, also
        write the active conda environment to a YML file in 'output_path'. Full
        simulation records are always written to the output '.pdb' or '.pdb.bz2'
        file(s), which can be used to reproduce any decoy without the scorefile.
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
    from dask.distributed import as_completed
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

from datetime import datetime
from pyrosetta.distributed.cluster.base import TaskBase, _get_residue_type_set
from pyrosetta.distributed.cluster.converters import (
    _parse_decoy_ids,
    _parse_environment,
    _parse_input_packed_pose,
    _parse_pyrosetta_build,
    _parse_scratch_dir,
    _parse_seeds,
    _parse_sha1,
    _parse_system_info,
    _parse_tasks,
)
from pyrosetta.distributed.cluster.initialization import _get_pyrosetta_init_args, _maybe_init_client
from pyrosetta.distributed.cluster.io import IO
from pyrosetta.distributed.cluster.logging_support import LoggingSupport
from pyrosetta.distributed.cluster.multiprocessing import user_spawn_thread
from pyrosetta.distributed.cluster.serialization import Serialization
from pyrosetta.distributed.cluster.utilities import SchedulerManager
from pyrosetta.distributed.cluster.validators import (
    _validate_dir,
    _validate_dirs,
    _validate_float,
    _validate_int,
    _validate_min_len,
)
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import (
    Any,
    List,
    NoReturn,
    Optional,
    TypeVar,
    Union,
)


G = TypeVar("G")


@attr.s(kw_only=True, slots=True, frozen=False)
class PyRosettaCluster(IO[G], LoggingSupport[G], SchedulerManager[G], TaskBase[G]):
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
    min_workers = attr.ib(
        type=int,
        default=1,
        validator=[_validate_int, attr.validators.instance_of(int)],
        converter=attr.converters.default_if_none(default=1),
    )
    max_workers = attr.ib(
        type=int,
        default=attr.Factory(
            lambda self: 1000 if (self.tasks_size < 1000) else self.tasks_size,
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
    scorefile_name = attr.ib(
        type=str,
        default="scores.json",
        validator=attr.validators.instance_of(str),
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
    environment = attr.ib(
        type=str,
        default=None,
        validator=attr.validators.instance_of(str),
        converter=_parse_environment,
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
                        "environment.yml",
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

    def __attrs_post_init__(self):
        _maybe_init_client()
        self._setup_logger()
        self._write_environment_file(self.environment_file)
        self.serializer = Serialization(compression=self.compression)
        self.clients_dict = self._setup_clients_dict()

    def distribute(self, *args: Any, protocols: Any = None, clients_indices: Any = None, resources: Any = None) -> Optional[NoReturn]:
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
                applied, the protocols will not run. See https://distributed.dask.org/en/stable/resources.html for more
                information.
                Default: None

        Returns:
            None
        """

        compressed_input_packed_pose = self.serializer.compress_packed_pose(self.input_packed_pose)
        protocols, protocol, seed, clients_index, resource = self._setup_protocols_protocol_seed(args, protocols, clients_indices, resources)
        clients, cluster, adaptive = self._setup_clients_cluster_adaptive()
        client_residue_type_set = _get_residue_type_set()
        extra_args = (
            self.decoy_ids,
            self.protocols_key,
            self.timeout,
            self.ignore_errors,
            self.logging_file,
            self.logging_level,
            self.DATETIME_FORMAT,
            self.compression,
            self.max_delay_time,
            client_residue_type_set,
        )
        seq = as_completed(
            [
                clients[clients_index].submit(
                    user_spawn_thread,
                    protocol,
                    compressed_input_packed_pose,
                    compressed_kwargs,
                    pyrosetta_init_kwargs,
                    *extra_args,
                    pure=False,
                    resources=resource,
                )
                for compressed_kwargs, pyrosetta_init_kwargs in (
                    self._setup_initial_kwargs(protocols, seed, task_kwargs) for task_kwargs in self.tasks
                )
                for _ in range(self.nstruct)
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
                kwargs = self.serializer.decompress_kwargs(compressed_kwargs)
                if not kwargs[self.protocols_key]:
                    self._save_results(compressed_packed_pose, kwargs)
                else:
                    if self.save_all:
                        self._save_results(
                            compressed_packed_pose,
                            self.serializer.deepcopy_kwargs(kwargs),
                        )
                    compressed_kwargs, pyrosetta_init_kwargs, protocol, clients_index, resource = self._setup_kwargs(
                        kwargs, clients_indices, resources
                    )
                    scatter = clients[clients_index].scatter(
                        (
                            protocol,
                            compressed_packed_pose,
                            compressed_kwargs,
                            pyrosetta_init_kwargs,
                            *extra_args,
                        ), hash=False,
                    )
                    seq.add(clients[clients_index].submit(user_spawn_thread, *scatter, pure=False, resources=resource))
                    self.tasks_size += 1
                    self._maybe_adapt(adaptive)

        self._maybe_teardown(clients, cluster)
        self._close_logger()


PyRosettaCluster.__doc__ = __doc__
