from .core import (
    register_builtin_container_traversal,
    dict_to_pose,
    dict_to_packed,
    pack_result,
    pose_result,
    to_packed,
    to_pose,
    to_dict,
    to_base64,
    to_pickle,
    PackedPose,
)

try:
    from .pandas import register_pandas_container_traversal

except ImportError:
    register_pandas_container_traversal = None


def register_container_traversal(generic_func, dict_func):
    register_builtin_container_traversal(generic_func, dict_func)
    if register_pandas_container_traversal:
        register_pandas_container_traversal(generic_func, dict_func)


register_container_traversal(to_pose, dict_to_pose)
register_container_traversal(to_packed, dict_to_packed)
register_container_traversal(to_dict, dict)
