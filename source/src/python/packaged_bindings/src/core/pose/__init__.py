__import__("rosetta.core.pose.datacache")

from __pose_all_at_once_ import *

from cStringIO import StringIO

from types import MethodType

import rosetta.utility
import rosetta.core.io.silent


def __pose_getstate__(self):
    """Generate a state dictionary containing a binary silent struct of the input pose."""

    silent_data = rosetta.core.io.silent.SilentFileData()
    silent_struct = rosetta.core.io.silent.BinaryProteinSilentStruct(self, "pickle")

    data = StringIO()
    silent_struct.print_header(rosetta.utility.ostream(data))
    silent_data.write_silent_struct(silent_struct, rosetta.utility.ostream(data))

    state = {
            "BinaryProteinSilentStruct" : data.getvalue(),
            "name" : self.pdb_info().name() if self.pdb_info() is not None else ""
            }

    return state

def __pose_setstate__(self, state):
    """Read state dictionary containing a binary silent struct into the pose."""

    data = state["BinaryProteinSilentStruct"]

    #Empty tag array reads all tags
    silent_data = rosetta.core.io.silent.SilentFileData()
    silent_data.read_stream(rosetta.utility.istream(StringIO(data)), rosetta.utility.vector1_string(), True)

    silent_struct = silent_data.get_structure("pickle")

    silent_struct.fill_pose(self)

    if state.get("name", ""):
        self.pdb_info().name(state["name"])

def __pose_reduce_ex__(self, protocol_level):
    """Pose reduce implementation."""
    return (Pose, tuple(), self.__getstate__())

# Add pickle support to pose objects.
Pose.__getstate__ = MethodType(__pose_getstate__, None, Pose)
Pose.__setstate__ = MethodType(__pose_setstate__, None, Pose)
Pose.__reduce_ex__ = MethodType(__pose_reduce_ex__, None, Pose)
