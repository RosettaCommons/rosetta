"""
The `PyRosettaCluster` Python class facilitates scalable and reproducible job distribution of user-defined
PyRosetta protocols efficiently parallelized on the user's local workstation, high-performance computing (HPC)
cluster, or elastic cloud computing infrastructure with available compute resources.

*Warning*: This class uses the `pickle` module to deserialize pickled `Pose` objects. Using the `pickle`
module is not secure, so please only run with input files you trust. Learn more about the `pickle` module
and its security `here <https://docs.python.org/3/library/pickle.html>`_.

Args:
    `tasks`:
        A `list` object of JSON-serializable `dict` objects, a callable returning an iterable of
        JSON-serializable `dict` objects, or a generator that yields JSON-serializable `dict` objects. During
        simulation execution, each task dictionary is automatically unpacked and collected by the variadic
        keyword parameter of each user-defined PyRosetta protocol; as a result, items may be accessed
        dynamically during PyRosetta protocol execution via the dictionary of keyword arguments, enabling
        dynamic control flow in PyRosetta protocols. In order to initialize PyRosetta with user-specified
        Rosetta command-line options at the start of each user-defined PyRosetta protocol, either or both of
        the `"extra_options"` and/or `"options"` keys may be defined in each task dictionary, where each value
        is either a `str`, `list`, or `dict` object representing the Rosetta command-line options.

        Default: `[{}]`

    `nstruct`:
        An `int` object specifying the number of repeats of the first user-defined PyRosetta protocol. The
        user can control the number of repeats of downstream PyRosetta protocols via returning multiple
        clones of any output decoys from any upstream PyRosetta protocols, or by cloning the input decoy
        multiple times inside a downstream PyRosetta protocol.

        Default: `1`

    `input_packed_pose`:
        An input `PackedPose` object that is accessible via the first positional-or-keyword parameter of the
        first user-defined PyRosetta protocol.

        Default: `None`

    `seeds`:
        A `list` of `int` objects specifying the PyRosetta pseudorandom number generator (RNG) seeds to use
        for each user-defined PyRosetta protocol. The length of the keyword argument value provided must be
        equal to the number of input user-defined PyRosetta protocols. Seeds are used in the same order that
        the user-defined PyRosetta protocols are executed.

        Default: `None`

    `decoy_ids`:
        A `list` of `int` objects specifying the decoy identification numbers to keep after executing the
        user-defined PyRosetta protocols. User-defined PyRosetta protocols may return an iterable of `Pose`
        and/or `PackedPose` objects, or yield `Pose` and/or `PackedPose` objects. To reproduce a particular
        decoy produced via the chain of user-provided PyRosetta protocols, the decoy number to keep for each
        protocol may be specified, where other decoys are discarded. Decoy numbers use zero-based indexing, so
        `0` is the first decoy generated from a particular user-defined PyRosetta protocol. The length of the
        keyword argument value provided must be equal to the number of input user-defined PyRosetta protocols,
        so that one decoy is kept for each user-defined PyRosetta protocol. Decoy identification numbers are
        applied in the same order that the user-defined PyRosetta protocols are executed.

        Default: `None`

    `client`:
        An initialized Dask `distributed.Client` object to be used as the Dask client interface to the local
        or remote Dask cluster. If `None`, then `PyRosettaCluster` initializes its own Dask client based on
        the `scheduler` keyword argument value. Deprecated by the `clients` keyword argument, but supported
        for legacy purposes. Either or both of the `client` or `clients` keyword argument values must be
        `None`.

        Default: `None`

    `clients`:
        A `list` or `tuple` object of initialized Dask `distributed.Client` objects to be used as the Dask
        client interface(s) to the local or remote Dask cluster(s). If `None`, then `PyRosettaCluster`
        initializes its own Dask client based on the `scheduler` keyword argument value. Optionally used in
        combination with the `PyRosettaCluster.distribute(clients_indices=...)` keyword argument. Either or
        both of the `client` or `clients` keyword argument values must be `None`. See the
        `PyRosettaCluster.distribute` method docstring for usage examples.

        Default: `None`

    `scheduler`:
        A `str` object of either `"sge"` or `"slurm"`, or `None`. If `"sge"`, then `PyRosettaCluster`
        schedules jobs using a `SGECluster` instance from the `dask-jobqueue` package. If `"slurm"`, then
        `PyRosettaCluster` schedules jobs using a `SLURMCluster` instance from the `dask-jobqueue` package.
        If `None`, then `PyRosettaCluster` schedules jobs using a `distributed.LocalCluster` instance. If
        `client` or `clients` keyword argument values are not `None`, then this keyword argument is ignored.

        Default: `None`

    `cores`:
        An `int` object specifying the total number of cores per job, which is passed to the
        `dask_jobqueue.SLURMCluster(cores=...)` or the `dask_jobqueue.SGECluster(cores=...)`
        keyword argument depending on the `scheduler` keyword argument value.

        Default: `1`

    `processes`:
        An `int` object specifying the total number of processes per job, which is passed to the
        `dask_jobqueue.SLURMCluster(processes=...)` or the `dask_jobqueue.SGECluster(processes=...)`
        keyword argument depending on the `scheduler` keyword argument value. This feature determines how
        many Python processes each Dask worker job will run.

        Default: `1`

    `memory`:
        A `str` object specifying the total amount of memory per job, which is input to the
        `dask_jobqueue.SLURMCluster(memory=...)` or the `dask_jobqueue.SGECluster(memory=...)` keyword
        argument depending on the `scheduler` keyword argument value.

        Default: `"4g"`

    `scratch_dir`:
        A `str` object specifying the absolute filesystem path to a scratch directory where temporary
        files may be written.

        Default: `"/temp"` if it exists, otherwise the current working directory.

    `min_workers`:
        An `int` object specifying the minimum number of workers to which to adapt during parallelization
        of user-defined PyRosetta protocols.

        Default: `1`

    `max_workers`:
        An `int` object specifying the maximum number of workers to which to adapt during parallelization
        of user-defined PyRosetta protocols.

        Default: `1000` if the number of user-defined task dictionaries passed to the `tasks` keyword argument
        value is `<1000`, otherwise the number of user-defined task dictionaries.

    `dashboard_address`:
        A `str` object specifying the port over which the Dask dashboard is forwarded. Particularly useful for
        diagnosing `PyRosettaCluster` performance in real-time.

        Default: `":8787"`

    `project_name`:
        A `str` object specifying the project name for this simulation. This keyword argument value is just
        added to the full simulation record for accounting purposes.

        Default: `datetime.now().strftime("%Y.%m.%d.%H.%M.%S.%f")` if not specified, else `"PyRosettaCluster"`
        if `None`.

    `simulation_name`:
        A `str` object specifying the particular name of this simulation. This keyword argument value is just
        added to the full simulation record for accounting purposes.

        Default: `project_name` if not specified, else `"PyRosettaCluster"` if `None`

    `output_path`:
        A `str` object specifying the absolute path of the output directory where the results will be written
        to disk. The directory will be created be created if it does not exist.

        Default: `"./outputs"`

    `output_decoy_types`:
        An iterable of `str` objects representing the output decoy filetypes to save during the simulation.
        Available options are: `".pdb"` for PDB files; `".pkl_pose"` for pickled `Pose` files; `".b64_pose"`
        for Base64-encoded pickled `Pose` files; and `".init"` for PyRosetta initialization files, each
        caching: the Rosetta command-line options (and PyRosetta initialization input files, if any)
        initialized on the head node, the `input_packed_pose` keyword argument value (if any), and an output
        decoy. Because each PyRosetta initialization file contains a copy of the PyRosetta initialization
        input files and input `PackedPose` object (if any), unless these objects are relatively small in size
        or there are relatively few expected output decoys, then it is recommended to run the
        `pyrosetta.distributed.cluster.export_init_file` function on only output decoys of interest after the
        simulation completes without specifying `".init"` in this iterable. If the `compressed` keyword
        argument value is set to `True`, then each output decoy file is further compressed by the `bzip2`
        library, and ".bz2" is automatically appended to the filename.

        Default: `[".pdb",]`

    `output_scorefile_types`:
        An iterable of `str` objects representing the output scorefile filetypes to save during the
        simulation. Available options are: `".json"` for a JSON Lines (JSONL)-formatted scorefile, and any
        filename extensions accepted by the `pandas.DataFrame.to_pickle(compression="infer")` method
        (including `".gz"`, `".bz2"`, and `".xz"`) for pickled `pandas.DataFrame` objects of scorefile data
        that may be analyzed using `pyrosetta.distributed.cluster.io.secure_read_pickle(compression="infer")`.
        Note that in order to write pickled `pandas.DataFrame` objects to disk, please ensure that
        `pyrosetta.secure_unpickle.add_secure_package("pandas")` has first been run. If using `pandas` version
        `>=3.0.0`, PyArrow-backed datatypes may be enabled by default; in this case, please ensure that
        `pyrosetta.secure_unpickle.add_secure_package("pyarrow")` has also first been run.

        See https://pandas.pydata.org/pdeps/0010-required-pyarrow-dependency.html and
        https://pandas.pydata.org/pdeps/0014-string-dtype.html for more information.

        Default: `[".json",]`

    `scorefile_name`:
        A `str` object specifying the name of the output JSON Lines (JSONL)-formatted scorefile, which must
        end in ".json". The scorefile location is always ```output_path` "/" `scorefile_name```. If `".json"`
        is not in the `output_scorefile_types` keyword argument value, the JSONL-formatted scorefile will not
        be output, but other scorefile types (if any) will still get the same filename stem (i.e., before the
        ".json" extension).

        Default: `"scores.json"`

    `simulation_records_in_scorefile`:
        A `bool` object specifying whether or not to write full simulation records to the output scorefile(s).
        If `True`, then write full simulation records to the output scorefile(s). This results in some
        redundant information in each entry, allowing downstream reproduction of a decoy of interest from the
        scorefile, but a larger scorefile storage footprint. If `False`, then write curtailed simulation
        records to the scorefile(s) containing only the `Pose.cache` dictionary data. This results in
        minimally redundant information in each entry, disallowing downstream reproduction of a decoy of
        interest from the output scorefile(s), but a smaller scorefile storage footprint. If `False`, also
        write the active Conda/Mamba/uv/Pixi environment configuration to a separate output file in the
        `output_path` keyword argument value for accounting purposes. Full simulation records are always
        written to the output decoy files, which may still be used to reproduce any decoy without the
        scorefile.

        Default: `False`

    `decoy_dir_name`:
        A `str` object specifying the directory name where the output decoy files will be saved. The directory
        location is always determined by: ```output_path` "/" `decoy_dir_name```.

        Default: `"decoys"`

    `logs_dir_name`:
        A `str` object specifying the directory name where the output log files will be saved. The directory
        location is always determined by: ```output_path` "/" `logs_dir_name```.

        Default: `"logs"`

    `logging_level`:
        A `str` object specifying the logging level of Python logging output to write to the log file.
        The available options are either: `"NOTSET"`, `"DEBUG"`, `"INFO"`, `"WARNING"`, `"ERROR"`, or
        `"CRITICAL"`. The output log file location is always determined by:
        ```output_path` "/" `logs_dir_name` "/" `simulation_name```.

        Default: `"INFO"`

    `logging_address`:
        A `str` object specifying the socket endpoint for sending and receiving log messages across a network,
        so log messages from user-defined PyRosetta protocols may be written to a single log file on the head
        node. The `str` object must take the format of a socket address (i.e.,
        "`<host address>`:`<port number>`") where the `<host address>` is either an IP address, "localhost", or
        Domain Name System (DNS)-accessible domain name, and the `<port number>` is a digit greater than or
        equal to "0". If the port number is "0", then the next free port number is selected.

        Default: `"localhost:0"` if the `scheduler` keyword argument value is `None` or either the `client` or
        `clients` keyword argument values specify instances of `distributed.LocalCluster`, else `"0.0.0.0:0"`.

    `compressed`:
        A `bool` object specifying whether or not to compress the output decoy files and output PyRosetta
        initialization files using the `bzip2` library, resulting in the appending of ".bz2" to output decoy
        files and PyRosetta initialization files. Also see the `output_decoy_types` and `output_init_file`
        keyword arguments.

        Default: `True`

    `compression`:
        A `str` object of either `"xz"`, `"zlib"` or `"bz2"`, or a `bool` or `None` object representing the
        internal compression library for pickled `Pose` objects and user-defined task dictionaries. The
        default of `True` uses `"xz"` for compression if it is installed, otherwise resorts to `"zlib"` for
        compression.

        Default: `True`

    `sha1`:
        A `str` or `None` object specifying the Git commit SHA-1 hash string of the local Git repository
        defining the simulation. If a non-empty `str` object is provided, then it is validated to match the
        Git commit SHA-1 hash string of the most recent commit in the local Git repository checked out in the
        current working directory, and then it is added to the simulation record for accounting. If an empty
        string is provided, then ensure that everything in the current working directory is committed to the
        local Git repository. If `None` is provided, then bypass SHA-1 hash string validation and set this
        value to an empty string.

        Default: `""`

    `ignore_errors`:
        A `bool` object specifying whether or not to ignore raised Python exceptions and thrown Rosetta
        segmentation faults in the user-defined PyRosetta protocols. This comes in handy when well-defined
        errors are sparse and sporadic (such as rare Rosetta segmentation faults), and the user would like
        `PyRosettaCluster` to continue running without otherwise raising a `WorkerError` exception.

        Default: `False`

    `timeout`:
        A `float` or `int` object specifying how many seconds to wait between `PyRosettaCluster` checking-in
        on the running user-defined PyRosetta protocols. If each PyRosetta protocol is expected to run quickly,
        then `0.1` seconds seems reasonable. If each PyRosetta protocol is expected to run slowly, then `>1`
        second seems reasonable.

        Default: `0.5`

    `max_delay_time`:
        A `float` or `int` object specifying the maximum number of seconds to sleep before returning the
        result(s) from each user-defined PyRosetta protocol back to the Dask client on the head node. If a
        Dask worker returns the result(s) from a PyRosetta protocol too quickly, the Dask scheduler needs to
        first register that the task is processing before it completes. In practice, in each PyRosetta
        protocol the runtime is subtracted from the `max_delay_time` keyword argument value, and the Dask
        worker sleeps for the remainder of the time (if any) before returning the result(s). It is recommended
        to set this option to at least 1 second, but longer times may be used as a safety throttle in cases
        of overwhelmed Dask scheduler processes. Because spawning a `billiard` subprocess for PyRosetta
        protocol execution may take ~3–5 seconds already before the PyRosetta protocol executes, this feature
        usually does not have an effect with the default value.

        Default: `3.0`

    `filter_results`:
        A `bool` object specifying whether or not to filter out empty `PackedPose` objects between
        user-defined PyRosetta protocols. When a PyRosetta protocol returns or yields `None`,
        `PyRosettaCluster` converts it to an empty `PackedPose` object that gets bound to the first
        positional-or-keyword parameter of the next PyRosetta protocol. If `True`, then filter out any empty
        `PackedPose` objects where there are no residues in the conformation as given by `PackedPose.empty()`.
        If `False`, then continue to pass the empty `PackedPose` objects to the next PyRosetta protocol. This
        is used for filtering out decoys mid-trajectory in-between PyRosetta protocols if PyRosetta protocols
        return or yield any `None`, empty `Pose`, or empty `PackedPose` objects.

        Default: `True`

    `save_all`:
        A `bool` object specifying whether or not to save all of the returned non-empty `PackedPose` objects
        from all user-defined PyRosetta protocols. This option may be used to checkpoint decoy trajectories
        after each PyRosetta protocol.

        Default: `False`

    `dry_run`:
        A `bool` object specifying whether or not to save output decoy files to disk. If `True`, then do not
        write output decoy files to disk. This feature may be useful for debugging.

        Default: `False`

    `norm_task_options`:
        A `bool` object specifying whether or not to normalize the task
        'options' and 'extra_options' values after PyRosetta initialization on the remote
        compute cluster. If `True`, then this enables more facile simulation reproduction
        by the use of the `ProtocolSettingsMetric` SimpleMetric to normalize the PyRosetta
        initialization options and by relativization of any input files and directory paths
        to the current working directory from which the task is running.

        Default: `True`

    `max_task_replicas`:
        An `int` or `None` object specifying the replication factor of tasks on Dask workers within the network
        (only via Dask's best effort). If an `int` object, the value must be greater than or equal to `0`. If
        `None`, then attempt to replicate all tasks on each Dask worker. Tasks are automatically deleted from
        each Dask worker upon task completion. Task replication improves resilience of the simulation when
        compute resources executing tasks are preempted midway through a user-defined PyRosetta protocol
        (e.g., due to using cloud spot instances or cluster backfill queues), so scattered data can be
        recovered. If a Dask worker is preempted during task execution, then the number of task retries is
        controlled by the Dask configuration parameter `distributed.scheduler.allowed-failures`, which may be
        manually configured prior to the simulation. Dask worker memory limits may also need to be increased to
        achieve the desired replication factor (see `memory` keyword argument). Using task replicas requires
        that either Dask's `ReduceReplicas` policy is disabled or that Dask's entire Active Memory Manager
        (AMM) is disabled, since replicated tasks consume additional memory per Dask worker. Task size in
        memory is dominated by the input `PackedPose` object; a rough estimate of additional memory usage is
        ~1 MB/task for a 500 residue protein. Task retries are only appropriate when PyRosetta protocols are
        side effect-free upon preemption, wherein tasks can be restarted without producing inconsistent
        external states if preempted midway through a PyRosetta protocol.

        See https://distributed.dask.org/en/stable/api.html#distributed.Client.replicate and
        https://docs.dask.org/en/stable/configuration.html for more information.

        Default: `0`

    `task_registry`:
        A `None` object or `str` object of either `"disk"` or `"memory"`. If `"disk"` is provided, then write
        the task registry to disk. If `"memory"` is provided, then keep the task registry in memory on the head
        node process. Maintaining a task registry improves the resilience of the simulation when compute
        resources executing tasks are preempted midway through a user-defined PyRosetta protocol (e.g., due to
        using cloud spot instances or cluster backfill queues); if scattered data cannot be recovered (see
        `max_task_replicas` keyword argument), then the task will be automatically resubmitted using the task
        input arguments cached in the task registry. If `"memory"` is provided, then task input arguments
        consume memory on the head node process, which is appropriate with fewer tasks (e.g., debugging
        pipelines). If `"disk"` is provided, then task input arguments consume disk space (in the `scratch_dir`
        keyword argument value), which is appropriate for production simulations. Task size is dominated by the
        input `PackedPose` object; a rough estimate of additional disk or memory usage is ~1 MB/task for a 500
        residue protein. Completed tasks are automatically deleted from the task registry upon task completion.
        If `None` is provided, then the task registry is not created, which is appropriate for non-preemptible
        compute resources. Task resubmissions are only appropriate when user-provided PyRosetta protocols are
        side effect-free upon preemption, wherein tasks can be restarted without producing inconsistent
        external states if preempted midway through a PyRosetta protocol.

        Default: `None`

    `cooldown_time`:
        A `float` or `int` object specifying how many seconds to sleep after the simulation is complete to
        allow loggers to flush. For very slow network filesystems, 2 or more seconds may be reasonable.

        Default: `0.5`

    `system_info`:
        A `dict` or `None` object specifying the system information and/or extra simulation informatio
        required to reproduce the simulation. If `None` is provided, then `PyRosettaCluster` automatically
        detects the platform and sets this value as the dictionary ``{"sys.platform": `sys.platform`}`` (e.g.,
        `{"sys.platform": "linux"}`). If a `dict` object is provided, then validate that the `"sys.platform"`
        key has a value equal to the current `sys.platform` result, and log a warning message if not.
        System information such as Amazon Machine Image (AMI) identifier and compute fleet instance type
        identifier may be stored in this dictionary, but it is not automatically validated upon reproduction
        simulations. This information is stored in the full simulation records for accounting.

        Default: `None`

    `pyrosetta_build`:
        A `str` or `None` object specifying the PyRosetta build signature as output by
        `pyrosetta._build_signature()`. If `None` is provided, then `PyRosettaCluster` automatically detects
        the PyRosetta build signature and sets this keyword argument value. If a non-empty `str` object is
        provided, then validate that the input PyRosetta build signature is equal to the active PyRosetta
        build signature, and raise an exception if not. This validation process ensures that reproduction
        simulations use an identical PyRosetta build signature from the original simulation. To bypass
        PyRosetta build signature validation with a warning message, an empty string ('') may be provided
        but does not assure reproducibility.

        Default: `None`

    `security`:
        A `bool` object or instance of `distributed.Security`, only having an effect if both `client=None` and
        `clients=None`, that is passed to Dask if using `scheduler=None` or passed to Dask-Jobqueue if using
        `scheduler="slurm"` or `scheduler="sge"`. If `True` is provided, then invoke the `cryptography` package
        to generate a `distributed.Security.temporary` object through Dask or Dask-Jobqueue. If a Dask
        `distributed.Security` object is provided, then pass it to Dask with `scheduler=None`, or pass it to
        Dask-Jobqueue with `scheduler="slurm"` or `scheduler="sge"` (where the `shared_temp_directory` keyword
        argument value of `SLURMCluster` or `SGECluster` is set to the `output_path` keyword argument value of
        `PyRosettaCluster`). If `False` is provided, then Dask TLS security is disabled regardless of the
        `scheduler` keyword argument value (which is *not* recommended for remote Dask clusters unless behind a
        trusted private network segment (i.e., a firewall). If `None` is provided, then `True` is used by
        default. In order to generate a `distributed.Security` object with the OpenSSL command-line interface,
        the `pyrosetta.distributed.cluster.generate_dask_tls_security` function may also be used (see docstring
        for more information) instead of the `cryptography` package.

        See https://distributed.dask.org/en/latest/tls.html#distributed.security.Security.temporary for more
        information.

        Default: `False` if `scheduler=None`, else `True`

    `max_nonce`:
        An `int` object greater than or equal to 1 specifying the maximum number of nonces to cache per process
        if Dask TLS security is disabled while using remote Dask clusters, which protects against replay
        attacks. If nonce caching is in use, each process (including the head node process and all Dask worker
        processes) cache nonces upon communication exchange over the network, which can increase memory usage
        in each process. A rough estimate of additional memory usage is ~0.2 KB per task per user-defined
        PyRosetta protocol per process. For example, submitting 1000 tasks with 2 PyRosetta protocols adds
        (~0.2 KB/task/protocol × 1000 tasks × 2 protocols) = ~0.4 MB of additional memory per processs. If
        memory usage per process permits, it is recommended to set this value to at least the number of tasks
        times the number of protocols submitted, so that every nonce from every communication exchange over the
        network gets cached.

        Default: `4096`

    `environment`:
        A `None` or `str` object specifying either the active Conda/Mamba environment YML file string, active
        uv project `uv.lock` file string, or active Pixi project `pixi.lock` file string. If `None` is
        provided, then generate an environment file string for the active Conda/Mamba/uv/Pixi environment and
        save it to the full simulation record. If a non-empty `str` object is provided, then validate it
        to match the active Conda/Mamba/uv/Pixi environment YML/lock file string and save it to the full
        simulation record. This ensures that reproduction simulations use an identical Conda/Mamba/uv/Pixi
        environment configuration to the original simulation. To bypass Conda/Mamba/uv/Pixi environment
        validation with a warning message, an empty string ('') may be provided, but does not assure
        reproducibility.

        Default: `None`

    `author`:
        A `str` object specifying the author(s) of the simulation that is written to the full simulation
        records and the output PyRosetta initialization file(s).

        Default: `""`

    `email`:
        A `str` object specifying the email address(es) of the author(s) of the simulation that is written to
        the full simulation records and the output PyRosetta initialization file(s).

        Default: `""`

    `license`:
        A `str` object specifying the license of the output data of the simulation that is written to the full
        simulation records and the output PyRosetta initialization file(s) (e.g., "ODC-ODbL", "CC BY-ND",
        "CDLA Permissive-2.0", etc.).

        Default: `""`

    `output_init_file`:
        A `str` object specifying the absolute path to the output PyRosetta initialization file that caches the
        `input_packed_pose` keyword argument value upon `PyRosettaCluster` instantiation. The file does not
        include any output decoys, and is optionally used for exporting PyRosetta initialization files with
        output decoys by the `pyrosetta.distributed.cluster.export_init_file` function after the simulation
        completes (see the `output_decoy_types` keyword argument). If `None` (or an empty `str` object (`""`))
        is provided, or `dry_run` keyword argument value is set to `True`, then skip writing an output ".init"
        file upon `PyRosettaCluster` instantiation. If skipped, it is recommended to run the
        `pyrosetta.dump_init_file` function before or after the simulation. If the `compressed` keyword
        argument value is set to `True`, then the output file is further compressed by the `bzip2` library,
        and ".bz2" is automatically appended to the filename.

        Default: ```output_path` "/" `project_name` "_" `simulation_name` "_pyrosetta.init"``

Returns:
    A `PyRosettaCluster` instance.
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
    import toolz
    from distributed import (
        Client,
        Future,
        Security,
        as_completed,
    )
    from distributed.scheduler import KilledWorker
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.core' requires the "
        + "third-party packages 'attrs', 'distributed', and 'toolz' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/attrs/\n"
        + "https://pypi.org/project/distributed/\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import logging
import os
import uuid

from concurrent.futures import CancelledError
from datetime import datetime
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import (
    Any,
    Dict,
    Generator,
    List,
    Optional,
    Tuple,
    TypeVar,
    Union,
)

from pyrosetta.distributed.cluster.base import (
    TaskBase,
    _get_residue_type_set,
)
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
from pyrosetta.distributed.cluster.exceptions import TaskCancelledError
from pyrosetta.distributed.cluster.hkdf import derive_task_key
from pyrosetta.distributed.cluster.initialization import (
    _get_pyrosetta_init_args,
    _maybe_init_client,
)
from pyrosetta.distributed.cluster.io import IO
from pyrosetta.distributed.cluster.logging_support import (
    LoggingSupport,
    MaskedBytes,
)
from pyrosetta.distributed.cluster.multiprocessing import user_spawn_thread
from pyrosetta.distributed.cluster.security import SecurityIO
from pyrosetta.distributed.cluster.serialization import (
    NonceCache,
    Serialization,
)
from pyrosetta.distributed.cluster.task_registry import (
    DiskTaskRegistry,
    MemoryTaskRegistry,
    UserArgs,
)
from pyrosetta.distributed.cluster.utilities import SchedulerManager
from pyrosetta.distributed.cluster.validators import (
    _validate_dir,
    _validate_dirs,
    _validate_float,
    _validate_int,
    _validate_logging_address,
    _validate_max_task_replicas,
    _validate_min_len,
    _validate_output_init_file,
    _validate_scorefile_name,
    _validate_tasks,
)

G = TypeVar("G")


@attr.s(kw_only=True, slots=True, frozen=False)
class PyRosettaCluster(IO[G], LoggingSupport[G], SchedulerManager[G], SecurityIO[G], TaskBase[G]):
    tasks = attr.ib(
        type=list,
        default=[{}],
        validator=[
            attr.validators.deep_iterable(
                member_validator=attr.validators.instance_of(dict),
                iterable_validator=attr.validators.instance_of(list),
            ),
            _validate_tasks,
        ],
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
        type=Optional[Client],
        default=None,
        validator=attr.validators.optional(
            attr.validators.instance_of(Client)
        ),
    )
    clients = attr.ib(
        type=Optional[List[Client]],
        default=None,
        validator=[
            attr.validators.optional(
                attr.validators.deep_iterable(
                    member_validator=attr.validators.instance_of(Client),
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
        validator=attr.validators.instance_of(str),
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
    max_task_replicas = attr.ib(
        type=Optional[int],
        default=0,
        validator=[
            _validate_max_task_replicas,
            attr.validators.optional(attr.validators.instance_of(int))
        ],
    )
    task_registry = attr.ib(
        type=Optional[str],
        default=None,
        validator=[
            attr.validators.optional(attr.validators.instance_of(str)),
            attr.validators.in_([None, "disk", "memory"]),
        ],
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
                            else "uv.lock"
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
    task_registry_dir = attr.ib(
        type=Optional[str],
        default=attr.Factory(
            lambda self: (
                os.path.join(self.scratch_dir, f"task_registry-{self.instance_id}")
                if self.task_registry == "disk"
                else None
            ),
            takes_self=True,
        ),
        init=False,
        validator=attr.validators.optional(attr.validators.instance_of(str)),
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
        """Pre-initialization hook for `PyRosettaCluster`."""
        _maybe_issue_environment_warnings()

    def __attrs_post_init__(self) -> None:
        """Post-initialization hook for `PyRosettaCluster`."""

        _maybe_init_client()
        self._setup_logger()
        self._cache_toml()
        self._write_environment_file(self.environment_file)
        self._write_init_file()
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
        if self.task_registry == "disk":
            self.registry = DiskTaskRegistry(
                instance_id=self.instance_id,
                compression=self.compression,
                task_registry_dir=self.task_registry_dir,
            )
        elif self.task_registry == "memory":
            self.registry = MemoryTaskRegistry(
                instance_id=self.instance_id,
                compression=self.compression,
            )
        else:
            self.registry = None

    def _get_submit_kwargs(
        self,
        resources: Optional[Dict[Any, Any]] = None,
        priority: Optional[int] = None,
        retries: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Setup `Client.submit` keyword arguments."""

        submit_kwargs = {"pure": False}
        # Omit resources keyword argument for distributed versions <2.1.0
        # or use default if user specifies `resources=None` in distributed versions >=2.1.0
        if resources is not None:
            submit_kwargs["resources"] = resources
        # Omit priority keyword argument for distributed versions <1.21.0
        # or use default if user specifies `priorities=None` in distributed versions >=1.21.0
        if priority is not None:
            submit_kwargs["priority"] = priority
        # Omit retries keyword argument for distributed versions <1.20.0
        # or use default if user specifies `retries=None` for distributed versions >=1.20.0
        if retries is not None:
            submit_kwargs["retries"] = retries

        return submit_kwargs

    def _create_future(
        self,
        client: Client,
        clients_index: int,
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
        """Scatter data and return submitted `user_spawn_thread` future."""

        task_id = uuid.uuid4().hex
        masked_key = MaskedBytes(derive_task_key(passkey, task_id))
        user_args = UserArgs(
            protocol_name=protocol_name,
            compressed_protocol=compressed_protocol,
            compressed_packed_pose=compressed_packed_pose,
            compressed_kwargs=compressed_kwargs,
            pyrosetta_init_kwargs=pyrosetta_init_kwargs,
            client_repr=repr(client),
            extra_args=extra_args,
            masked_key=masked_key,
            task_id=task_id,
        )
        scatter = client.scatter(user_args, broadcast=False, hash=False)
        if self.max_task_replicas != 0:
            client.replicate(scatter, n=self.max_task_replicas)
        submit_kwargs = self._get_submit_kwargs(resources=resource, priority=priority, retries=retry)
        future = client.submit(user_spawn_thread, scatter, **submit_kwargs)
        if self.task_registry:
            self.registry.set(
                future.key,
                clients_index=clients_index,
                user_args=user_args,
                submit_kwargs=submit_kwargs,
            )

        return future

    def _recreate_future(
        self,
        client: Client,
        clients_index: int,
        user_args: UserArgs,
        submit_kwargs: Dict[str, Any],
    ) -> Future:
        """Re-scatter data and return submitted 'user_spawn_thread' future."""

        scatter = client.scatter(user_args, broadcast=False, hash=False)
        if self.max_task_replicas != 0:
            client.replicate(scatter, n=self.max_task_replicas)
        future = client.submit(user_spawn_thread, scatter, **submit_kwargs)
        if self.task_registry:
            self.registry.set(
                future.key,
                clients_index=clients_index,
                user_args=user_args,
                submit_kwargs=submit_kwargs,
            )

        return future

    def _run(
        self,
        *args: Any,
        protocols: Any = None,
        clients_indices: Any = None,
        resources: Any = None,
        priorities: Any = None,
        retries: Any = None,
    ) -> Generator[Tuple[PackedPose, Dict[Any, Any]], None, None]:
        """
        Run user-provided PyRosetta protocols on a local or remote compute cluster using
        the user-customized PyRosettaCluster instance. Either arguments or the 'protocols'
        keyword argument is required. If both are provided, then the 'protocols' keyword
        argument value gets concatenated after the input arguments.

        *Warning*: This method uses the `cloudpickle` and `pickle` modules to serialize and deserialize `Pose`
        objects, arbitrary Python types in `Pose.cache` dictionaries, `pandas.DataFrame` objects (if
        configured), user-defined task dictionaries, user-defined PyRosetta protocols, and other user-provided
        data. Using the `cloudpickle` and `pickle` modules is not secure, so please only run this method with
        input data you fully understand and trust. Learn more about the `cloudpickle` and `pickle` modules and
        their security `here <https://github.com/cloudpipe/cloudpickle>`_ and
        `here <https://docs.python.org/3/library/pickle.html>`_.

        Examples:

        Basic usage:
            >>> PyRosettaCluster().distribute(protocol_1)
            >>> PyRosettaCluster().distribute(protocols=protocol_1)
            >>> PyRosettaCluster().distribute(protocol_1, protocol_2, protocol_3)
            >>> PyRosettaCluster().distribute(protocols=(protocol_1, protocol_2, protocol_3))
            >>> PyRosettaCluster().distribute(protocol_1, protocol_2, protocols=[protocol_3, protocol_4])

        Run with two Dask clients:
            >>> # Run `protocol_1` on `client_1`,
            >>> # then `protocol_2` on `client_2`,
            >>> # then `protocol_3` on `client_1`,
            >>> # then `protocol_4` on `client_2`:
            >>> PyRosettaCluster(clients=[client_1, client_2]).distribute(
            ...     protocols=[protocol_1, protocol_2, protocol_3, protocol_4],
            ...     clients_indices=[0, 1, 0, 1],
            ... )

        Run with multiple Dask clients:
            >>> # Run `protocol_1` on `client_2`,
            >>> # then `protocol_2` on `client_3`,
            >>> # then `protocol_3` on `client_1`:
            >>> PyRosettaCluster(clients=[client_1, client_2, client_3]).distribute(
            ...     protocols=[protocol_1, protocol_2, protocol_3],
            ...     clients_indices=[1, 2, 0],
            ... )

        Run with one Dask client and compute resource constraints:
            >>> # Run `protocol_1` on `client_1` with Dask worker resource constraints "GPU=2",
            >>> # then `protocol_2` on `client_1` with Dask worker resource constraints "MEMORY=100e9",
            >>> # then `protocol_3` on `client_1` without Dask worker resource constraints:
            >>> PyRosettaCluster(client=client_1).distribute(
            ...     protocols=[protocol_1, protocol_2, protocol_3],
            ...     resources=[{"GPU": 2}, {"MEMORY": 100e9}, None],
            ... )

        Run with two Dask clients and compute resource constraints:
            >>> # Run `protocol_1` on `client_1` with Dask worker resource constraints "GPU=2",
            >>> # then `protocol_2` on `client_2` with Dask worker resource constraints "MEMORY=100e9":
            >>> PyRosettaCluster(clients=[client_1, client_2]).distribute(
            ...     protocols=[protocol_1, protocol_2],
            ...     clients_indices=[0, 1],
            ...     resources=[{"GPU": 2}, {"MEMORY": 100e9}],
            ... )

        Run with task priorities:
            >>> # Run protocols with depth-first task execution:
            >>> PyRosettaCluster().distribute(
            ...     protocols=[protocol_1, protocol_2, protocol_3, protocol_4],
            ...     priorities=[0, 10, 20, 30],
            ... )

        Run with task retries:
            >>> # Run protocols with up to three retries per failed task during `protocol_3` and `protocol_4`:
            >>> PyRosettaCluster(ignore_errors=False).distribute(
            ...     protocols=[protocol_1, protocol_2, protocol_3, protocol_4],
            ...     retries=[0, 0, 3, 3],
            ... )

        Args:
            `*args`:
                Optional callables of type `types.GeneratorType` or `types.FunctionType` representing
                user-defined PyRosetta protocols in the order to be executed.

            `protocols`:
                An iterable of extra callable user-defined PyRosetta protocols; i.e., an iterable of objects
                of `types.GeneratorType` and/or `types.FunctionType` types, or a single callable of type
                `types.GeneratorType` or `types.FunctionType`.

                Default: `None`

            `clients_indices`:
                A `list` or `tuple` object of `int` objects, where each `int` object represents a zero-based
                index corresponding to the initialized Dask `distributed.Client` object(s) passed to the
                `PyRosettaCluster(clients=...)` keyword argument value. If not `None`, then the length of the
                `clients_indices` object must equal the number of protocols passed to the
                `PyRosettaCluster.distribute` method.

                Default: `None`

            `resources`:
                A `list` or `tuple` object of `dict` objects, where each `dict` object represents an abstract,
                arbitrary resource to constrain which Dask workers execute the user-defined PyRosetta protocols.
                If `None`, then do not impose resource constaints on any PyRosetta protocols. If not `None`,
                then the length of the `resources` object must equal the number of PyRosetta protocols passed
                to the `PyRosettaCluster.distribute` method, such that each resource specified indicates the
                unique resource constraints for the protocol at the corresponding index of the PyRosetta
                protocols passed to the `PyRosettaCluster.distribute` method. Note that this feature is only
                useful when one passes in their own instantiated Dask client(s) with Dask workers set up with
                various resource constraints. If Dask workers were not instantiated to satisfy the specified
                resource constraints, PyRosetta protocols will hang indefinitely by design because the Dask
                scheduler is waiting for Dask workers that meet the specified resource constraints so that it
                may schedule these tasks. Unless Dask workers were created with these resource tags applied,
                the PyRosetta protocols will not run.

                See https://distributed.dask.org/en/stable/resources.html for more information.

                Default: `None`

            `priorities`:
                A `list` or `tuple` object of `int` objects, where each `int` object sets the Dask scheduler
                priority for the corresponding user-defined PyRosetta protocol (i.e., indexed the same as the
                `client_indices` keyword argument value). If `None`, then no explicit priorities are set. If
                not `None`, then the length of this value must equal the number of PyRosetta protocols passed
                to the `PyRosettaCluster.distribute` method, and each `int` value determines the Dask scheduler
                priority for the tasks applied to that PyRosetta protocol.

                Breadth-first task execution (default):
                    When all user-defined PyRosetta protocols have an identical priority (e.g.,
                    `[0] * len(protocols)` or `None`), then all tasks enter the Dask scheduler's queue with
                    equal priority. Under equal priority, Dask mainly schedules tasks in a first-in, first-out
                    behavior. When Dask worker resources are saturated, this causes all tasks submitted to
                    upstream PyRosetta protocols to run to completion before tasks are scheduled to execute
                    downstream PyRosetta protocols, producing a breadth-first task execution behavior across
                    PyRosetta protocols.

                Depth-first task execution:
                    To allow tasks to run through all user-defined PyRosetta protocols before all tasks
                    applied to upstream PyRosetta protocols complete, assign increasing priorities to
                    downstream protocols (e.g., `list(range(0, len(protocols) * 10, 10))`). Once a task
                    completes an upstream PyRosetta protocol, it is applied to the next downstream PyRosetta
                    protocol with a higher priority than tasks still queued for upstream PyRosetta protocols,
                    so tasks may run through all user-defined PyRosetta protocols to completion as Dask worker
                    resources become available. This produces a depth-first task execution behavior across
                    PyRosetta protocols when Dask worker resources are saturated.

                See https://distributed.dask.org/en/stable/priority.html for more information.

                Default: `None`

            `retries`:
                A `list` or `tuple` of `int` objects, where each `int` object (≥0) sets the number of allowed
                automatic retries of each failed task that was applied to the corresponding user-defined
                PyRosetta protocol (i.e., indexed the same as `client_indices` keyword argument value). If an
                `int` object (≥0) is provided, then apply that number of allowed automatic retries to all
                PyRosetta protocols. If `None` is provided, then no explicit retries are allowed. If not `None`
                and not an `int` object, then the length of this value must equal the number of PyRosetta
                protocols passed to the `PyRosettaCluster.distribute` method, and each `int` value determines
                the number of automatic retries the Dask scheduler allows for that the tasks applied to that
                PyRosetta protocol. Allowing retries of failed tasks may be useful if the PyRosetta protocol
                raises a standard Python exception or Rosetta throws a segmentation fault in the `billiard`
                subprocess while the Dask worker remains alive and `PyRosettaCluster(ignore_errors=False)` is
                configured. If `PyRosettaCluster(ignore_errors=True)` is configured, then protocols failing due
                to standard Python exceptions or Rosetta segmentation faults will still be considered
                successes, and this keyword argument has no effect since these PyRosetta protocol errors are
                ignored. Note that if a compute resource executing a PyRosetta protocol is preempted, then the
                Dask worker process does not remain alive and the Dask scheduler registers that failed task as
                incomplete or cancelled. In this case, the number of allowed task retries is controlled by the
                Dask configuration parameter `distributed.scheduler.allowed-failures`; please use the
                `max_task_replices` and `task_registry` keyword arguments of `PyRosettaCluster` for further
                configuration of task retries after compute resource preemption.

                See https://distributed.dask.org/en/latest/scheduling-state.html#task-state for more
                information.

                Default: `None`
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
                    clients_index,
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
            except CancelledError as ex:
                if not self.task_registry:
                    raise TaskCancelledError(
                        future.key, "Please enable the task registry to resubmit this task."
                    ) from ex
                _task_record_values = self.registry.get(future.key)
                if _task_record_values is None:
                    raise TaskCancelledError(
                        future.key, "Task arguments could not be recovered from the task registry."
                    ) from ex
                logging.info(
                    f"Caught exception {type(ex).__name__}: {ex}. Resubmitting task from task registry and continuing."
                )
                clients_index, user_args, submit_kwargs = _task_record_values
                seq.add(
                    self._recreate_future(
                        clients[clients_index],
                        clients_index,
                        user_args,
                        submit_kwargs,
                    )
                )
                self.tasks_size += 1
                self._maybe_adapt(adaptive)
                continue
            except KilledWorker as ex:
                logging.error(
                    f"Caught exception {type(ex).__name__}: {ex}. Dropping task and continuing."
                )
                continue
            finally:
                if self.task_registry:
                    self.registry.pop(future.key)
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
                            clients_index,
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

        if self.task_registry:
            self.registry.clear()
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
    ) -> Generator[Tuple[PackedPose, Dict[str, Any]], None, None]:
        # See `generate.__doc__` explicitly set below

        if self.sha1 != "":
            logging.warning(
                "Use of the `PyRosettaCluster.generate` method for reproducible simulations is not supported! "
                + "`PyRosettaCluster` reproduces decoys from the output files written to disk. Subsequent code run "
                + "on these results is not being saved by `PyRosettaCluster`. Use the `PyRosettaCluster.distribute` "
                + "method for reproducible simulations. To silence this warning and continue without using version "
                + "control, set the `sha1` keyword argument value of `PyRosettaCluster` to `None`."
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
    ) -> None:
        # See `distribute.__doc__` explicitly set below

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
        remote compute cluster and optionally post-processing or visualizing output `PackedPose` objects in
        memory. Importantly, subsequent code run on the yielded results is not captured by `PyRosettaCluster`,
        and so use of this method does not ensure reproducibility of the simulation. Use the
        `PyRosettaCluster.distribute` method for reproducible simulations.

        Each yielded result is a `tuple` object with a `PackedPose` object as the first element and a `dict`
        object as the second element. The `PackedPose` object represents a returned or yielded `PackedPose`
        (or `Pose` or `NoneType`) object from the most recently executed user-defined PyRosetta protocol. The
        `dict` object represents the optionally returned or yielded dictionary of keyword arguments from the
        same most recently executed PyRosetta protocol (see the `protocols` keyword argument). If
        `PyRosettaCluster(save_all=True)` is used, tuples are yielded after each PyRosetta protocol, otherwise
        tuples are yielded after the final PyRosetta protocol. Tuples are yielded in the order in which they
        arrive back to the Dask client(s) from the distributed cluster (which may differ from the order that
        tasks are submitted, due to tasks running asynchronously). If `PyRosettaCluster(dry_run=True)` is used,
        then tuples are still yielded but output decoy files are not written to disk.

        See https://docs.dask.org/en/latest/futures.html#distributed.as_completed for more information.

        Extra examples:

        Iterate over results in real-time as they are yielded from the cluster:
            >>> for packed_pose, kwargs in PyRosettaCluster().generate(protocols):
            ...     ...

        Iterate over submissions to the same Dask client:
            >>> client = Client()
            >>> for packed_pose, kwargs in PyRosettaCluster(client=client).generate(protocols):
            ...     # Post-process results on head node asynchronously from results generation
            ...     prc = PyRosettaCluster(
            ...         input_packed_pose=packed_pose,
            ...         client=client,
            ...         logs_dir_name=f"logs_{uuid.uuid4().hex}", # Make sure to write new log files
            ...     )
            ...     for packed_pose, kwargs in prc.generate(other_protocols):
            ...         ...

        Iterate over two PyRosettaCluster instances, each managing one Dask client, creating additional
        overhead:
            >>> client_1 = Client()
            >>> client_2 = Client()
            >>> for packed_pose, kwargs in PyRosettaCluster(client=client_1).generate(protocols):
            ...     # Post-process results on head node asynchronously from results generation
            ...     prc = PyRosettaCluster(
            ...         input_packed_pose=packed_pose,
            ...         client=client_2,
            ...         logs_dir_name=f"logs_{uuid.uuid4().hex}", # Make sure to write new log files
            ...     )
            ...     for packed_pose, kwargs in prc.generate(other_protocols):
            ...         ...

        Iterate over one PyRosettaCluster instance managing two Dask clients, reducing overhead:
            >>> # Using multiple `distributed.as_completed` iterators on the head node creates additional
            >>> # overhead. If post-processing on the head node is not required between user-defined PyRosetta
            >>> # protocols, the preferred method is to distribute PyRosetta protocols within a single
            >>> # `PyRosettaCluster.generate` method call using the `clients_indices` keyword argument:
            >>> prc_generate = PyRosettaCluster(clients=[client_1, client_2]).generate(
            ...     protocols=[protocol_1, protocol_2],
            ...     clients_indices=[0, 1],
            ...     resources=[{"GPU": 1}, {"CPU": 1}],
            ... )
            ... for packed_pose, kwargs in prc_generate:
            ...     # Post-process results on head node asynchronously from results generation
            ...     ...

        Yields:
            ``(`PackedPose`, `dict`)`` tuples from the most recently executed user-defined PyRosetta protocol
            if `PyRosettaCluster(save_all=True)` is used, otherwise from the final PyRosetta protocol.
        """

PyRosettaCluster.__doc__ = __doc__
