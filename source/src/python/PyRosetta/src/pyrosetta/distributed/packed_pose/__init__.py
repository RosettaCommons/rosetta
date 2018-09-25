from .core import (
    register_builtin_container_traversal,
    dict_to_pose,
    dict_to_packed,
    pack_result,
    to_packed,
    to_pose,
    to_dict,
    PackedPose,
)

register_builtin_container_traversal(to_pose, dict_to_pose)
register_builtin_container_traversal(to_packed, dict_to_packed)
register_builtin_container_traversal(to_dict, dict)

try:
    from .pandas import register_pandas_container_traversal

    register_pandas_container_traversal(to_pose, dict_to_pose)
    register_pandas_container_traversal(to_packed, dict_to_packed)
    register_pandas_container_traversal(to_dict, dict)

except ImportError:
    pass
