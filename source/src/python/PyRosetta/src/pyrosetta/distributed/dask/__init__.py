# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

import pyrosetta
import pyrosetta.distributed


def _normflags(flags):
    """Normalize tuple/list/str of flags into str."""
    if not isinstance(flags, str):
        flags = " ".join(flags)
    return " ".join(" ".join([line.split("#")[0] for line in flags.split("\n")]).split())


def worker_extra(init_flags=None, local_directory=None):
    """Format flags and local directory for dask worker preload.

    Convenience function that optionally accepts unformatted PyRosetta command
    line flags as the `init_flags` option, and a working directory for dask
    worker processes to output files as the `local_directory` option, which are
    then formatted into a dask worker process extra commands list. If
    `init_flags` is not specified, PyRosetta is initialized on the dask worker
    process without command line flags. If `local_directory` is not specified,
    the dask worker default directory ./dask-worker-space is used.

    For example, to execute on a SLURM scheduler via dask_jobqueue:

    from dask_jobqueue import SLURMCluster
    import pyrosetta.distributed.dask
    cluster = SLURMCluster(
        cores=1,
        processes=1,
        memory="4GB",
        extra=pyrosetta.distributed.dask.worker_extra(init_flags, local_directory)
    )
    """

    extras = []
    if local_directory:
        extras.extend(["--local-directory", local_directory])

    if init_flags:

        # Specify flags with preceeding ' ' so that command line parser passes
        # concatenated list of flags into worker preload handler.
        extras.extend(
            [
                "--preload pyrosetta.distributed.dask.worker ' {0}'".format(
                    _normflags(init_flags)
                )
            ]
        )

    return extras


def init_notebook(init_flags=None):
    """Initialize PyRosetta with command line flags.

    For example, to execute in a Jupyter Notebook:

    import pyrosetta.distributed.dask
    pyrosetta.distributed.dask.init_notebook(init_flags)
    """
    kwargs = {}

    if init_flags:
        kwargs["extra_options"] = _normflags(init_flags)
    else:
        kwargs["extra_options"] = ""

    pyrosetta.distributed.maybe_init(**kwargs)
