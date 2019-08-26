# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.


import functools
import logging
import threading
import random
import traceback

import pyrosetta
import pyrosetta.rosetta.basic.random

_logger = logging.getLogger("pyrosetta.distributed")

__all__ = ["init", "maybe_init", "requires_init", "with_lock"]

# Access lock for any non-threadsafe calls into rosetta internals.
# Intended to provide a threadsafe api surface area to `distributed`.
_access_lock = threading.RLock()


def _normflags(flags):
    """Normalize tuple/list/str of flags into str."""
    if not isinstance(flags, str):
        flags = " ".join(flags)
    return " ".join(" ".join([line.split("#")[0] for line in flags.split("\n")]).split())
    

def with_lock(func):
    """Function decorator that protects access to rosetta internals."""
    @functools.wraps(func)
    def fwrap(*args, **kwargs):
        # Use DEBUG - 1 to work around errors with logging level init
        # under nosetests
        if _logger.isEnabledFor(logging.DEBUG-1):
            tb = traceback.extract_stack()[-5:-2]
            _logger.log(logging.DEBUG-1, "entering with_lock:\n%s", "\n".join(traceback.format_list(tb)))
        _logger.debug("with_lock call: %s", func)
        with _access_lock:
            try:
                return func(*args, **kwargs)
            finally:
                _logger.debug("with_lock finished: %s", func)

    return fwrap


@with_lock
def maybe_init(**kwargs):
    """Call pyrosetta.init if not already initialized."""
    # if "set_logging_handler" not in kwargs:
        # kwargs["set_logging_handler"] = "logging"
    #TODO In multithreaded builds the python-based logging handler appears to fail...
    if "extra_options" not in kwargs:
        kwargs["extra_options"] = "-out:levels all:warning"
    if "silent" not in kwargs:
        kwargs["silent"] = True

    if not pyrosetta.rosetta.basic.was_init_called():
        _logger.info("maybe_init performing pyrosetta initialization: %s", kwargs)
        pyrosetta.init(**kwargs)

    # Thread init culled from multithreaded job distributer
    # Maybe converge into a "thread_init" function in rosetta core?
    rgen = pyrosetta.rosetta.numeric.random.rg()
    if not rgen.initialized():
        _logger.info("maybe_init performing per-thread rng initialization")
        options = pyrosetta.rosetta.basic.options.process()
        rgs = pyrosetta.rosetta.basic.random.RandomGeneratorSettings()
        rgs.initialize_from_options(options)
        pyrosetta.rosetta.basic.random.init_random_generators(random.randint(0, 0xFFFFFF), rgs.rng_type())


def requires_init(func):
    """Asserting that pyrosetta is initialized before function execution."""

    @functools.wraps(func)
    def fwrap(*args, **kwargs):
        maybe_init()

        return func(*args, **kwargs)

    return fwrap


def init(options=None, **kwargs):
    """Initialize PyRosetta with command line options."""
    if options and ("extra_options" not in kwargs):
        kwargs["extra_options"] = _normflags(options)
    maybe_init(**kwargs)
