# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"

import pyrosetta
import pyrosetta.distributed


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
    ```
    from dask_jobqueue import SLURMCluster
    import pyrosetta.distributed.dask
    cluster = SLURMCluster(
        cores=1,
        processes=1,
        memory="4GB",
        worker_extra_args=pyrosetta.distributed.dask.worker_extra(init_flags, local_directory)
    )
    ```
    Set `extra` instead of `worker_extra_args` if using dask-jobqueue version < 0.8.0
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
                    pyrosetta.distributed._normflags(init_flags)
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
    pyrosetta.distributed.init(options=init_flags)
