# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

import click
import pyrosetta.distributed


@click.command()
@click.argument("init_flags", nargs=1, type=click.STRING)
def dask_setup(worker, init_flags=()):
    """dask-worker '--preload' Command Line Interface Creation Kit (click)
    hook to initialize dask-worker with PyRosetta command line flags.
    """
    flags = " ".join("".join(init_flags).split())
    kwargs = {"extra_options": flags}
    pyrosetta.distributed.maybe_init(**kwargs)
