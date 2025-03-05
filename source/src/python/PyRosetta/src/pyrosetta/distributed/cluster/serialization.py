# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import attr
import bz2
import logging
import pyrosetta.distributed.io as io
import sys
import zlib

try:
    import lzma as xz
except ImportError:
    pass

try:
    import cloudpickle
    import toolz
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.serialization' requires the "
        + "third-party packages 'cloudpickle' and 'toolz' as dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/cloudpickle/\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

from functools import singledispatch, wraps
from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import (
    Any,
    Dict,
    Callable,
    Generic,
    NoReturn,
    Optional,
    TypeVar,
    Union,
)

T = TypeVar("T", bound=Callable[..., Any])
G = TypeVar("G")


def _parse_compression(obj: Any) -> Optional[Union[str, bool]]:
    """Parse the input `compression` attribute of `Serialization` class."""
    _error_msg = (
        "The 'compression' keyword argument parameter must be one of the `str` objects "
        + "'xz', 'zlib', or 'bz2', or an object of type `bool` or `NoneType`. Received: '{0}'"
    )

    @singledispatch
    def converter(obj: Any) -> NoReturn:
        raise ValueError(_error_msg.format(obj))

    @converter.register(str)
    def _is_str(obj: str) -> str:
        if obj not in ("xz", "zlib", "bz2"):
            raise ValueError(_error_msg.format(obj))

        return obj

    @converter.register(bool)
    def _is_bool(obj: bool) -> Union[bool, str]:
        if obj == True:
            if "lzma" in sys.modules:
                _compression = "xz"
            else:
                logging.warning("Could not import 'xz' library... Using 'zlib' for compression.")
                _compression = "zlib"
            return _compression
        else:
            return obj

    @converter.register(type(None))
    def _is_none(obj: None) -> str:
        return obj

    return converter(obj)


def update_scores(packed_pose: PackedPose) -> PackedPose:
    """
    Cache scores into the `PackedPose` object that are not cached in the `Pose` object and do not
    have keys with reserved scoretypes, then return the updated `PackedPose` object.

    Args:
        packed_pose: the input `PackedPose` object in which to update scores.

    Returns:
        A new `PackedPose` object with scores cached in its `Pose` object if scores could be cached,
        otherwise the input `PackedPose` object.
    """
    _pose = packed_pose.pose
    _pose_scoretypes = set(_pose.scores.keys())
    _reserved_scoretypes = _pose.__scores_accessor._reserved.union(_pose_scoretypes)
    _filtered_scores = toolz.dicttoolz.keyfilter(
        lambda scoretype: scoretype not in _reserved_scoretypes,
        packed_pose.scores,
    )
    if _filtered_scores:
        packed_pose = packed_pose.update_scores(_filtered_scores)

    return packed_pose


@attr.s(kw_only=False, slots=False, frozen=False)
class Serialization(Generic[G]):
    compression = attr.ib(
        type=Optional[Union[str, bool]],
        default="xz",
        validator=attr.validators.optional(attr.validators.instance_of((str, bool))),
        converter=_parse_compression,
    )

    def __attrs_post_init__(self):
        if self.compression == "xz":
            if "lzma" not in sys.modules:
                raise ImportError(
                    (
                        "Using 'xz' for compression requires installing the 'xz' package into your python environment. "
                        + "For installation instructions, visit:\n"
                        + "https://anaconda.org/anaconda/xz\n"
                    )
                )
            self.encoder = xz.compress
            self.decoder = xz.decompress
        elif self.compression == "zlib":
            self.encoder = self.zlib_compress
            self.decoder = zlib.decompress
        elif self.compression == "bz2":
            self.encoder = bz2.compress
            self.decoder = bz2.decompress
        elif self.compression in (False, None):
            self.encoder = None
            self.decoder = None

    @classmethod
    def zlib_compress(cls, obj):
        return zlib.compress(obj, 9)

    def requires_compression(func):
        @wraps(func)
        def wrapper(self, obj):
            if all(x is not None for x in (self.encoder, self.decoder)):
                return func(self, obj)
            else:
                logging.debug("Compression/decompression is disabled.")
                return obj

        return wrapper

    @requires_compression
    def compress_packed_pose(self, packed_pose: Any) -> Union[NoReturn, None, bytes]:
        """
        Compress a `PackedPose` object with the `bz2` module. If the 'packed_pose' argument parameter
        is `None`, then just return `None`.

        Args:
            packed_pose: the input `PackedPose` object to compress. If `None`, then just return `None`.

        Returns:
            A `bytes` object representing the compressed `PackedPose` object, or a `NoneType` object.

        Raises:
            `TypeError` if the 'packed_pose' argument parameter is not of type `NoneType` or `PackedPose`.
        """
        if packed_pose is None:
            compressed_packed_pose = None
        elif isinstance(packed_pose, PackedPose):
            packed_pose = update_scores(packed_pose)
            compressed_packed_pose = self.encoder(packed_pose.pickled_pose)
        else:
            raise TypeError(
                "The 'packed_pose' argument parameter must be of type `NoneType` or `PackedPose`."
            )

        return compressed_packed_pose

    @requires_compression
    def decompress_packed_pose(self, compressed_packed_pose: Any) -> Union[NoReturn, None, PackedPose]:
        """
        Decompress a `bytes` object with the `bz2` and `cloudpickle` modules. If the 'compressed_packed_pose'
        argument parameter is `None`, then just return `None`.

        Args:
            compressed_packed_pose: the input `bytes` object to decompress. If `None`, then just return `None`.

        Returns:
            A `PackedPose` object representing the decompressed `bytes` object, or a `NoneType` object.

        Raises:
            `TypeError` if the 'compressed_packed_pose' argument parameter is not of type `NoneType` or `bytes`.
        """
        if compressed_packed_pose is None:
            packed_pose = None
        elif isinstance(compressed_packed_pose, bytes):
            pose = cloudpickle.loads(self.decoder(compressed_packed_pose))
            packed_pose = io.to_packed(pose)
        else:
            raise TypeError(
                "The 'compressed_packed_pose' argument parameter must be of type `NoneType` or `bytes`."
            )

        return packed_pose

    @requires_compression
    def compress_kwargs(self, kwargs: Any) -> Union[NoReturn, bytes]:
        """
        Compress a `dict` object with the `cloudpickle` and `bz2` modules.

        Args:
            kwargs: the input `dict` object to compress.

        Returns:
            A `bytes` object representing the compressed `dict` object.

        Raises:
            `TypeError` if the 'kwargs' argument parameter is not of type `dict`.
        """
        if isinstance(kwargs, dict):
            return self.encoder(cloudpickle.dumps(kwargs))
        else:
            raise TypeError("The 'kwargs' argument parameter must be of type `dict`.")

    @requires_compression
    def decompress_kwargs(self, compressed_kwargs: Any) -> Union[NoReturn, Dict[Any, Any]]:
        """
        Decompress a `bytes` object with the `bz2` and `cloudpickle` modules.

        Args:
            compressed_kwargs: the input `bytes` object to decompress.

        Returns:
            A `dict` object representing the decompressed `bytes` object.

        Raises:
            `TypeError` if the 'compressed_packed_pose' argument parameter is not of type `bytes`.
        """
        if isinstance(compressed_kwargs, bytes):
            return cloudpickle.loads(self.decoder(compressed_kwargs))
        else:
            raise TypeError("The 'compressed_kwargs' argument parameter must be of type `bytes`.")

    @classmethod
    def deepcopy_kwargs(cls, kwargs: Any) -> Union[NoReturn, Dict[Any, Any]]:
        """
        The `cloudpickle` module makes it possible to serialize Python constructs not supported
        by the default `pickle` module from the Python standard library.

        Args:
            kwargs: the `dict` object to be deep copied.

        Returns:
            A deep copy of the `dict` object.

        Raises:
            `TypeError` if the 'kwargs' argument parameter is not of type `dict`.
        """
        if isinstance(kwargs, dict):
            return cloudpickle.loads(cloudpickle.dumps(kwargs))
        else:
            raise TypeError("The 'kwargs' argument parameter must be of type `dict`.")
