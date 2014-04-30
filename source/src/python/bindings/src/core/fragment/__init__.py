from _core_fragment_ import *

from cStringIO import StringIO

from types import MethodType

import rosetta.utility

def __fragment_getstate__(self):
    """Generate a state dictionary containing a serialized fragment frame containing the source fragment."""

    fragment_frame = Frame()
    fragment_frame.add_fragment(self)

    data = StringIO()
    fragment_frame.show(rosetta.utility.ostream(data))

    return {"FragmentFrame" : data.getvalue()}

def __fragment_setstate__(self, state):
    """Read state dictionary containing a serialized fragment frame containing the source fragment."""

    data = state["FragmentFrame"]

    frame_list = FrameList()

    fragment_io = FragmentIO()
    fragment_io.read_data( rosetta.utility.istream(StringIO(data)), frame_list )

    assert frame_list.flat_size() == 1

    self.copy(frame_list.fragID(1).fragment())

def __fragment_reduce_ex__(self, protocol_level):
    """Fragment reduce implementation, generates annotated fragments.""" 
    return (AnnotatedFragData, ("pickle", 0), self.__getstate__())

FragData.__getstate__ = MethodType(__fragment_getstate__, None, FragData)
FragData.__setstate__ = MethodType(__fragment_setstate__, None, FragData)
FragData.__reduce_ex__ = MethodType(__fragment_reduce_ex__, None, FragData)

AnnotatedFragData.__getstate__ = MethodType(__fragment_getstate__, None, AnnotatedFragData)
AnnotatedFragData.__setstate__ = MethodType(__fragment_setstate__, None, AnnotatedFragData)
AnnotatedFragData.__reduce_ex__ = MethodType(__fragment_reduce_ex__, None, AnnotatedFragData)
