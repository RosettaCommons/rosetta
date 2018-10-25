import warnings

import logging
logger = logging.getLogger(__name__)

try:
    import blosc

    # Activate blosc GIL-release mode, allowing fork-space multitheaded decompression
    # https://github.com/Blosc/bloscpack/issues/58
    # https://github.com/Blosc/c-blosc/pull/241
    logger.debug("blosc.set_releasegil(True)")
    blosc.set_releasegil(True)
except ImportError:
    logger.debug("ImportError loading blosc, falling back to uncompressed pickle archives.")
    blosc = None
    pass

try:
    import pyrosetta.rosetta.cereal as cereal
except ImportError:
    logger.debug("ImportError loading pyrosetta.rosetta.cereal, " 
                 "non-serialization build will not support pickling.")
    cereal = None
    pass

import pyrosetta.rosetta as rosetta
import pyrosetta.rosetta.utility as utility

from pyrosetta.distributed import with_lock

@with_lock
def __cereal_getstate__(self):
    # Defer blosc import to getstate to allow *optional* dependency
    if not cereal:
        raise NotImplementedError(
            "__cereal_getstate__ requires pyrosetta '--serialization' build.")

    oss = rosetta.std.ostringstream()
    self.save(cereal.BinaryOutputArchive(oss))

    result = {
        "cereal_archive_version": utility.Version.version()
    }

    if blosc:
        result["blosc_cereal_binary_archive"] = blosc.compress(oss.bytes())
    else:
        warnings.warn(
            "__cereal_getstate__ without 'blosc' installation. "
            "Archive size will be prohibitively large. "
            "Install (pip) `blosc` and/or (conda) `python-blosc`."
        )
        result["cereal_binary_archive"] = oss.bytes()

    return result

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
