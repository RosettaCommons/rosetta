import click
import pyrosetta.distributed


@click.command()
@click.argument('init_flags', nargs=1, type=click.STRING)
def dask_setup(worker, init_flags=()):
    """Dask worker --preload hook to init with PyRosetta command line flags."""
    flags = " ".join("".join(init_flags).split())
    kwargs = {"extra_options": flags}
    pyrosetta.distributed.maybe_init(**kwargs)
