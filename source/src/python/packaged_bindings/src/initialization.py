import sys
import os
import os.path

def rosetta_database_from_env():
    """Read rosetta database directory from environment or standard install locations.

    Returns database path if found, else None."""

    from rosetta import logger
    from rosetta import __path__ as rosetta_root

    candidate = os.path.join(rosetta_root[0], "database")
    if os.path.isdir(candidate):
        database = os.path.abspath(candidate)
        logger.info('Found rosetta database at: %s', database)
        return database

    # No database found.
    return None

def init(options='-ex1 -ex2aro', extra_options='', set_logging_handler=True):
    """Initialize Rosetta.  Includes core data and global options.

    options string with default Rosetta command-line options args.
            (default: '-ex1 -ex2aro')
    kargs -
        extra_options - Extra command line options to pass rosetta init.
                        (default None)
        set_logging_handler - Route rosetta tracing through logging logger 'rosetta'.
                        (default True)

    Examples:
        init()                     # uses default flags
        init(extra_options='-pH')  # adds flags to supplement the default
        init('-pH -database /home/me/pyrosetta/rosetta_database')  # overrides default flags - be sure to include the dB last
    """

    # Initialize rosetta logging handlers.
    from rosetta import logger, __version__

    import rosetta.logging_support as logging_support

    logging_support.initialize_logging()

    if set_logging_handler:
        logging_support.set_logging_handler()

    # Setup library-level exit callbacks
    from rosetta.utility import PythonPyExitCallback
    PythonPyExitCallback.setup()

    # Perform core/protocol level initialization
    from rosetta.utility import vector1_string

    args = ['pyrosetta'] + options.split() + extra_options.split()

    # Attempt to resolve database location from environment if not present, else fallback
    # to rosetta's standard resolution
    if not "-database" in args:
        database = rosetta_database_from_env()
        if database is not None:
            args.extend(["-database", database])

    init_args = vector1_string()
    init_args.extend(args)

    logger.info("Version: %s", __version__)

    try:
        from rosetta.protocols.init import init
    except ImportError:
        from rosetta.core.init import init

    init(init_args)
