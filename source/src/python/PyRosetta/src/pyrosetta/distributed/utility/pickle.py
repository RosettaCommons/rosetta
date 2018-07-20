import logging

try:
    import blosc
except ImportError:
    blosc = None
    pass

try:
    import pyrosetta.rosetta.cereal as cereal
except ImportError:
    cereal = None
    pass

import pyrosetta.rosetta as rosetta
import pyrosetta.rosetta.utility as utility

from pyrosetta.distributed import with_lock

logger = logging.getLogger(__name__)

@with_lock
def __cereal_getstate__(self):
    # Defer blosc import to getstate to allow *optional* dependency
    if not cereal:
        raise NotImplementedError(
            "__cereal_getstate__ requires pyrosetta '--serialization' build.")
    if not blosc:
        raise NotImplementedError(
            "__cereal_getstate__ requires 'blosc' installation.")

    oss = rosetta.std.ostringstream()
    self.save(cereal.BinaryOutputArchive(oss))
    return {
        "blosc_cereal_binary_archive": blosc.compress(oss.bytes()),
        "cereal_archive_version": utility.Version.version()
    }

@with_lock
def __cereal_setstate__(self, state):
    if not cereal:
        raise NotImplementedError(
            "__cereal_setstate__ requires pyrosetta '--serialization' build.")

    self.__init__()
    try:
        if "blosc_cereal_binary_archive" in state and blosc:
            iss = rosetta.std.istringstream(blosc.decompress(state["blosc_cereal_binary_archive"]))
        elif "cereal_binary_archive" in state:
            iss = rosetta.std.istringstream(state["cereal_binary_archive"])
        else:
            if "blosc_cereal_binary_archive" in state:
                raise ValueError(
                    "No blosc, unable to load compressed pickle state: %s" % tuple(state.keys()))
            else:
                raise ValueError(
                    "Unable to load unknown pickle state: %s" % tuple(state.keys()))

        self.load(cereal.BinaryInputArchive(iss))
    except Exception:
        logger.exception(
            "Error unpickling ceral archive type: %r"
            " archive_version: %r current_version: %r",
            type(self),
            state.get("cereal_archive_version", None),
            utility.Version.version())
        raise

def set_cereal_pickleable(klass):
    klass.__getstate__ = __cereal_getstate__
    klass.__setstate__ = __cereal_setstate__
