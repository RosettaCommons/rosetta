import sys

if sys.version_info.major >= 3:
    from functools import singledispatch
    import collections.abc
else:
    from pkgutil import simplegeneric as singledispatch

import pickle
import base64

import pyrosetta.rosetta.core.pose as pose
import pyrosetta.distributed

__all__ = ["pack_result", "to_packed", "to_pose", "to_dict", "PackedPose"]


class PackedPose:
    """Serializable, read-only access to a serialized pose object.

    PackedPose functions as the equivalent of the "Pose" object within the pyrosetta.distributed
    namespace. It holds a serialized pose object and a cached dictionary of scores extracted
    from the serialized pose. It can be "unpacked" into a working pose or isomorphicly represented
    as primitive datatypes or pandas data structures.

    The "primitive" datatype represention of PackedPosed is a dict containing, at least, the
    key "packed_pose" holding a base64-encoded unicode string of the pickled pose object. This
    representation supports transparent interconversion between json and the packed format. In
    contrast, PackedPose includes the *bytes* representation of the pose object and is not suitable
    for serialization as text.

    It should be noted that all pickled representations are *highly* compressible.
    """

    __slots__ = ("scores", "pickled_pose")

    def __init__(self, pose_or_pack):
        """Create a packed pose from pose, pack, or pickled bytes."""
        if isinstance(pose_or_pack, pose.Pose):
            self.pickled_pose = pickle.dumps(pose_or_pack)
            self.scores = dict(pose_or_pack.scores)

        elif isinstance(pose_or_pack, PackedPose):
            self.pickled_pose = pose_or_pack.pickled_pose
            self.scores = pose_or_pack.scores

        elif isinstance(pose_or_pack, bytes):
            self.pickled_pose = pose_or_pack
            self.scores = {}

        else:
            raise ValueError("Unknown input type.", type(pose_or_pack))

    @property
    @pyrosetta.distributed.requires_init
    def pose(self):
        return pickle.loads(self.pickled_pose)


def pack_result(func):
    @functools.wraps(func)
    def wrap(*args, **kwargs):
        return to_packed(func(*args, **kwargs))

    return wrap


@singledispatch
def to_packed(pose_or_pack):
    return PackedPose(pose_or_pack)


@to_packed.register(type(None))
def none_to_packed(none):
    return None


@to_packed.register(str)
def str_to_packed(b64_encoded_pickle):
    return to_packed(base64.b64decode(b64_encoded_pickle, validate=True))


def dict_to_packed(packed_dict):
    pack = to_packed(packed_dict["pickled_pose"])

    pack.scores = dict(packed_dict)
    pack.scores.pop("pickled_pose")

    return pack


@singledispatch
def to_pose(pack):
    return PackedPose(pack).pose


@to_pose.register(type(None))
def none_to_pose(none):
    return None


@to_pose.register(pose.Pose)
def pose_to_pose(p):
    return p


@to_pose.register(str)
def str_to_pose(b64_encoded_pickle):
    return to_pose(base64.b64decode(b64_encoded_pickle, validate=True))


def dict_to_pose(packed_dict):
    return to_pose(packed_dict["pickled_pose"])


@singledispatch
def to_dict(pose_or_pack):
    pack = PackedPose(pose_or_pack)

    result = {}
    result.update(pack.scores)
    result["pickled_pose"] = base64.b64encode(pack.pickled_pose).decode()

    return result


@to_dict.register(type(None))
def none_to_dict(none):
    return None


def register_builtin_container_traversal(generic_func, dict_func):
    @generic_func.register(dict)
    def dict_traversal(maybe_packed_dict):
        if "pickled_pose" in maybe_packed_dict:
            return dict_func(maybe_packed_dict)
        else:
            return {k: generic_func(v) for k, v in maybe_packed_dict.items()}

    @generic_func.register(list)
    @generic_func.register(tuple)
    @generic_func.register(set)
    def container_traveral(container):
        return container.__class__(map(generic_func, container))
        
    if sys.version_info.major >= 3:
        @generic_func.register(collections.abc.Generator)
        def generator_traversal(generator):
            return (generic_func(v) for v in generator)

    return generic_func
