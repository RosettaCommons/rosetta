# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    import attr
    import cloudpickle
    import msgpack
    import toolz
    from distributed import get_worker
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.serialization' requires the "
        + "third-party packages 'attrs', 'cloudpickle', 'distributed', 'msgpack', "
        + "and 'toolz' as dependencies!\nPlease install these packages into your "
        + "python environment. For installation instructions, visit:\n"
        + "https://pypi.org/project/attrs/\n"
        + "https://pypi.org/project/cloudpickle/\n"
        + "https://pypi.org/project/distributed/\n"
        + "https://pypi.org/project/msgpack/\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import bz2
import logging
import os
import pyrosetta.distributed.io as io
import sys
import warnings
import zlib

try:
    import lzma as xz
except ImportError:
    pass

from collections import deque
from functools import partial, singledispatch, wraps
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.secure_unpickle import SecureSerializerBase
from typing import (
    Any,
    Dict,
    Callable,
    Generic,
    NoReturn,
    Optional,
    TypeVar,
    Union,
    cast,
)

from pyrosetta.distributed.cluster.hkdf import MaskedBytes, compare_digest, hmac_digest


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
    _pose_scoretypes = set(_pose.cache.keys())
    _reserved_scoretypes = _pose.__cache_accessor._reserved.union(_pose_scoretypes)
    _filtered_scores = toolz.dicttoolz.keyfilter(
        lambda scoretype: scoretype not in _reserved_scoretypes,
        packed_pose.scores,
    )
    if _filtered_scores:
        warnings.warn(
            "PyRosettaCluster detected scores in the `PackedPose.scores` dictionary that have "
            + "not been added to the `Pose.cache` dictionary! PyRosettaCluster is now adding "
            + f"the following scores to the `Pose.cache` dictionary: {tuple(_filtered_scores.keys())}. "
            + "Please use the `packed_pose = packed_pose.update_scores(**scores_dict)` syntax instead of "
            + "the `packed_pose.scores[key] = value` syntax to suppress this warning.",
            SyntaxWarning,
            stacklevel=2,
        )
        packed_pose = packed_pose.update_scores(_filtered_scores)

    return packed_pose


@attr.s(kw_only=False, slots=True, frozen=True)
class MessagePacking(Generic[G]):
    """PyRosettaCluster MessagePack base class."""
    pack = attr.ib(
        type=partial,
        default=partial(msgpack.packb, use_bin_type=True),
        init=False,
        validator=attr.validators.instance_of(partial),
    )
    unpack = attr.ib(
        type=partial,
        default=partial(msgpack.unpackb, raw=False),
        init=False,
        validator=attr.validators.instance_of(partial),
    )


@attr.s(kw_only=True, slots=False, frozen=False)
class NonceCache(Generic[G]):
    """PyRosettaCluster nonce cache base class."""
    instance_id = attr.ib(
        type=str,
        validator=attr.validators.instance_of(str),
    )
    prk = attr.ib(
        type=Optional[bytes],
        default=None,
        repr=False,
        converter=lambda k: MaskedBytes(k) if isinstance(k, bytes) else k,
        validator=attr.validators.optional(attr.validators.instance_of(MaskedBytes)),
    )
    max_nonce = attr.ib(
        type=int,
        validator=attr.validators.instance_of(int),
    )
    _seen = attr.ib(
        type=set,
        default=attr.Factory(set, takes_self=False),
        init=False,
        validator=attr.validators.instance_of(set),
    )
    _order = attr.ib(
        type=deque,
        default=attr.Factory(
            lambda self: deque(maxlen=self.max_nonce),
            takes_self=True,
        ),
        init=False,
        validator=attr.validators.instance_of(deque),
    )
    _debug = attr.ib(
        type=bool,
        default=False,
        init=False,
        validator=attr.validators.instance_of(bool),
    )

    def __attrs_pre_init__(self) -> None:
        MessagePacking.__init__(self)

    @staticmethod
    def _get_state(self) -> Dict[str, Any]:
        """
        A method used to override the default `NonceCache.__getstate__()` method
        that sets the pseudo-random key value to `None` in the returned state.
        """
        state = self.__dict__.copy()
        state["prk"] = None

        return state

    def _cache_nonce(self, sealed: bytes) -> None:
        """
        Run Hash-based Message Authentication Code (HMAC) verification and cache nonces
        for replay protection without data decompression.
        """
        package = self.unpack(sealed)
        if not isinstance(package, dict) or package.get("v", None) != 1:
            _err_msg = (
                "Invalid sealed package or version on {0} nonce cache! "
                + f"Received: {type(package)!r}"
            )
            try:
                _worker = get_worker()
                raise ValueError(_err_msg.format("worker"))
            except BaseException:
                raise ValueError(_err_msg.format("host"))

        _instance_id = package["a"] # `str`: PyRosettaCluster instance identifier/App
        _data = package["d"] # `bytes`: Data bytestring
        _mac = package["m"] # `bytes`: MAC (HMAC tag)
        _nonce = package["n"] # `bytes`: Nonce
        _version = package["v"] # `int`: Version

        if _instance_id != self.instance_id:
            raise ValueError("PyRosettaCluster instance identifier mismatch in sealed package.")

        msg = self.pack([_instance_id, _data, _nonce, _version])
        _expected_mac = hmac_digest(bytes(self.prk), msg)
        if not compare_digest(_expected_mac, _mac):
            _err_msg = (
                "Task HMAC verification failed during nonce cache on {0}!\n"
                + f"Expected: {_expected_mac!r}\n"
                + f"Value:    {_mac!r}\n"
            )
            try:
                _worker = get_worker()
                raise SystemExit(_err_msg.format("worker"))
            except BaseException:
                raise SystemExit(_err_msg.format("host"))

        if _nonce is not None:
            if _nonce in self._seen:
                # Replay protection
                _err_msg = (
                    "PyRosettaCluster detected a repeat nonce on the {0} for the instance identifier "
                    + f"'{self.instance_id}', which might indicate a replay attack is in progress! "
                    + "Exiting process for security. Please ensure that `PyRosettaCluster(security=True)` "
                    + f"is enabled in future PyRosettaCluster simulations. Received: '{_nonce}'."
                )
                try:
                    _worker = get_worker()
                    raise SystemExit(_err_msg.format("worker"))
                except BaseException:
                    raise SystemExit(_err_msg.format("host"))
            self._seen.add(_nonce)
            self._order.append(_nonce)
            while len(self._seen) > self._order.maxlen:
                _expired = self._order.popleft()
                self._seen.discard(_expired)

            if self._debug:
                _memory_usage = round(sum(map(sys.getsizeof, (self._seen, self._order))) / 1e3, 3)  # KB
                _msg = (
                    f"Size={len(self._seen)}/{self.max_nonce}; "
                    + f"Memory usage: {_memory_usage} KB; "
                    + f"Example: {sorted(self._seen)[0]}"
                )
                try:
                    _worker = get_worker()
                    print(f"Remote worker ({_worker.contact_address}) nonce cache: {_msg}")
                except:
                    print(f"Local host nonce cache: {_msg}")


@attr.s(kw_only=True, slots=False, frozen=False)
class Serialization(Generic[G]):
    """PyRosettaCluster serialization base class."""
    instance_id = attr.ib(
        type=Optional[str],
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(str)),
    )
    prk = attr.ib(
        type=Optional[bytes],
        default=None,
        repr=False,
        converter=lambda k: MaskedBytes(k) if isinstance(k, bytes) else k,
        validator=attr.validators.optional(attr.validators.instance_of(MaskedBytes)),
    )
    compression = attr.ib(
        type=Optional[Union[str, bool]],
        default="xz",
        validator=attr.validators.optional(attr.validators.instance_of((str, bool))),
        converter=_parse_compression,
    )
    with_nonce = attr.ib(
        type=bool,
        default=False,
        validator=attr.validators.instance_of(bool),
    )
    _nonce_size = attr.ib(
        type=int,
        default=32, # bytes
        init=False,
        validator=attr.validators.instance_of(int),
    )


    def __attrs_pre_init__(self) -> None:
        MessagePacking.__init__(self)

    def __attrs_post_init__(self) -> None:
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

    def __getstate__(self) -> Dict[str, Any]:
        """
        A method used to override the default `Serialization.__getstate__()` method
        that sets the pseudo-random key value to `None` in the returned state.
        """
        state = self.__dict__.copy()
        state["prk"] = None

        return state

    @classmethod
    def zlib_compress(cls, obj: bytes) -> bytes:
        """Compress an object with `zlib` level 9."""
        return zlib.compress(obj, 9)

    def requires_compression(func: T) -> T:
        """
        Wrapper testing if compression is enabled, and skips compression if it's disabled.
        """
        @wraps(func)
        def wrapper(self, obj: Any) -> Any:
            if all(x is not None for x in (self.encoder, self.decoder)):
                return func(self, obj)
            else:
                logging.debug("Compression/decompression is disabled.")
                return obj

        return cast(T, wrapper)

    def _seal(self, data: bytes) -> bytes:
        """Seal data with MessagePack."""
        if self.instance_id is None or self.prk is None:
            raise ValueError(
                "Sealing requires the 'instance_id' and 'prk' instance attributes."
            )

        version = 1
        nonce = os.urandom(self._nonce_size) if self.with_nonce else None
        msg = self.pack([self.instance_id, data, nonce, version])
        mac = hmac_digest(self.prk, msg)
        package = {
            "a": self.instance_id, # `str`: PyRosettaCluster instance identifier/App
            "d": data, # `bytes`: Data bytestring
            "m": mac, # `bytes`: MAC (HMAC tag)
            "n": nonce, # `bytes`: Nonce
            "v": version, # `int`: Version
        }

        return self.pack(package)

    def _unseal(self, obj: bytes) -> bytes:
        """
        Unseal data with MessagePack and perform Hash-based
        Message Authentication Code (HMAC) verification.
        """
        if self.instance_id is None or self.prk is None:
            raise ValueError(
                "Unsealing requires the 'instance_id' and 'prk' instance attributes."
            )

        package = self.unpack(obj)
        if not isinstance(package, dict) or package.get("v", None) != 1:
            raise ValueError(f"Invalid sealed package or version. Received: {type(package)}")

        _instance_id = package["a"] # `str`: PyRosettaCluster instance identifier/App
        _data = package["d"] # `bytes`: Data bytestring
        _mac = package["m"] # `bytes`: MAC (HMAC tag)
        _nonce = package["n"] # `bytes`: Nonce
        _version = package["v"] # `int`: Version

        if self.instance_id is not None and _instance_id != self.instance_id:
            raise ValueError(
                "Instance identifier mismatch in packaged data!\n"
                + f"'{_instance_id}' != '{self.instance_id}'"
            )

        msg = self.pack([_instance_id, _data, _nonce, _version])
        _expected_mac = hmac_digest(self.prk, msg)
        if not compare_digest(_expected_mac, _mac):
            raise SystemExit(
                "Task HMAC verification failed upon unsealing the data!\n"
                + f"Expected: {_expected_mac!r}\n"
                + f"Value:    {_mac!r}\n"
            )

        return _data

    @requires_compression
    def compress_packed_pose(self, packed_pose: Any) -> Union[NoReturn, None, bytes]:
        """
        Compress a `PackedPose` object with the custom serialization module. If the 'packed_pose' argument
        parameter is `None`, then just return `None`.

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
            compressed_packed_pose = self.encoder(io.to_pickle(packed_pose))
        else:
            raise TypeError(
                "The 'packed_pose' argument parameter must be of type `NoneType` or `PackedPose`."
            )

        return compressed_packed_pose

    @requires_compression
    def decompress_packed_pose(self, compressed_packed_pose: Any) -> Union[NoReturn, None, PackedPose]:
        """
        Decompress a `bytes` object with the custom serialization module and secure implementation of the
        `pickle` module. If the 'compressed_packed_pose' argument parameter is `None`, then just return `None`.

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
            packed_pose = io.to_packed(SecureSerializerBase.secure_loads(self.decoder(compressed_packed_pose)))
        else:
            raise TypeError(
                "The 'compressed_packed_pose' argument parameter must be of type `NoneType` or `bytes`."
            )

        return packed_pose

    def loads_object(self, compressed_obj: bytes) -> Any:
        """Unseal data and run the `cloudpickle.loads` method."""
        data = self._unseal(compressed_obj)
        buffer = self.decoder(data) if self.decoder else data

        return cloudpickle.loads(buffer)

    def dumps_object(self, obj: Any) -> bytes:
        """Run the `cloudpickle.dumps` method and seal data."""
        pickled = cloudpickle.dumps(obj)
        buffer = self.encoder(pickled) if self.encoder else pickled

        return self._seal(buffer)

    def compress_kwargs(self, kwargs: Any) -> Union[NoReturn, bytes]:
        """
        Compress a `dict` object with the `cloudpickle` and custom serialization modules.

        Args:
            kwargs: the input `dict` object to compress.

        Returns:
            A `bytes` object representing the compressed `dict` object.

        Raises:
            `TypeError` if the 'kwargs' argument parameter is not of type `dict`.
        """
        if isinstance(kwargs, dict):
            return self.dumps_object(kwargs)
        else:
            raise TypeError("The 'kwargs' argument parameter must be of type `dict`.")

    def decompress_kwargs(self, compressed_kwargs: bytes) -> Union[NoReturn, Dict[Any, Any]]:
        """
        Decompress a `bytes` object with the custom serialization and `cloudpickle` modules.

        Args:
            compressed_kwargs: the input `bytes` object to decompress.

        Returns:
            A `dict` object representing the decompressed `bytes` object.

        Raises:
            `TypeError` if the 'compressed_packed_pose' argument parameter is not of type `bytes`.
            `TypeError` if the returned kwargs is not of type `dict`.
        """
        if isinstance(compressed_kwargs, bytes):
            kwargs = self.loads_object(compressed_kwargs)
            if not isinstance(kwargs, dict):
                raise TypeError(f"Decoded kwargs must be of type `dict`. Received: {type(kwargs)}")
            return kwargs
        else:
            raise TypeError("The 'compressed_kwargs' argument parameter must be of type `bytes`.")

    def compress_object(self, obj: Any) -> bytes:
        """
        Compress an object with the `cloudpickle` and custom serialization modules.

        Args:
            obj: the input object to compress.

        Returns:
            A `bytes` object representing the compressed object.
        """
        return self.dumps_object(obj)

    def decompress_object(self, compressed_obj: bytes) -> Any:
        """
        Decompress a `bytes` object with the custom serialization and `cloudpickle` modules.

        Args:
            compressed_obj: the input `bytes` object to decompress.

        Returns:
            An object representing the decompressed `bytes` object.

        Raises:
            `TypeError` if the 'compressed_obj' argument parameter is not of type `bytes`.
        """
        if isinstance(compressed_obj, bytes):
            return self.loads_object(compressed_obj)
        else:
            raise TypeError("The 'compressed_obj' argument parameter must be of type `bytes`.")

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
