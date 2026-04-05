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
        + "virtual environment. For installation instructions, visit:\n"
        + "https://pypi.org/project/attrs/\n"
        + "https://pypi.org/project/cloudpickle/\n"
        + "https://pypi.org/project/distributed/\n"
        + "https://pypi.org/project/msgpack/\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

try:
    import lzma as xz
except ImportError:
    pass

import bz2
import logging
import os
import pyrosetta.distributed.io as io
import sys
import warnings
import zlib

from collections import deque
from functools import (
    partial,
    singledispatch,
    wraps,
)
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.secure_unpickle import SecureSerializerBase
from typing import (
    AbstractSet,
    Any,
    Deque,
    Dict,
    Callable,
    NoReturn,
    Optional,
    TypeVar,
    Union,
    cast,
)

from pyrosetta.distributed.cluster.hkdf import (
    MaskedBytes,
    compare_digest,
    hmac_digest,
)
from pyrosetta.distributed.cluster.type_defs import CallableType


def _parse_compression(obj: Any) -> Optional[Union[str, bool]]:
    """Parse the `compression` keyword argument value of the `Serialization` class."""

    _error_msg = (
        "The `compression` keyword argument value must be one of the `str` objects "
        + "'xz', 'zlib', or 'bz2', or an object of type `bool` or `None`. Received: '{0}'"
    )

    @singledispatch
    def converter(obj: Any) -> NoReturn:
        raise ValueError(_error_msg.format(type(obj)))

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
    def _is_none(obj: None) -> None:
        return obj

    return converter(obj)


def update_scores(packed_pose: PackedPose) -> PackedPose:
    """
    Cache scores from the `PackedPose.scores` dictionary that are not cached in the `Pose.cache` dictionary and
    do not have keys with reserved scoretypes, then return a new `PackedPose` object.

    *Warning*: This function uses the `pickle` module to deserialize pickled `Pose` objects. Using the `pickle`
    module is not secure, so please only run with input files you trust. Learn more about the `pickle` module
    and its security `here <https://docs.python.org/3/library/pickle.html>`_.

    Args:
        `packed_pose`:
            An input `PackedPose` object in which to update scores.

    Returns:
        A new `PackedPose` object, with scores cached in its `Pose.cache` dictionary if scores could be cached.
    """

    _pose = packed_pose.pose
    _pose_scoretypes = set(_pose.cache.all_keys)
    _reserved_scoretypes = _pose.cache._reserved.union(_pose_scoretypes)
    _filtered_scores = toolz.dicttoolz.keyfilter(
        lambda scoretype: scoretype not in _reserved_scoretypes,
        packed_pose.scores,
    )
    if _filtered_scores:
        warnings.warn(
            "`PyRosettaCluster` detected scores in the `PackedPose.scores` dictionary that have "
            + "not been added to the `Pose.cache` dictionary! `PyRosettaCluster` is now adding "
            + f"the following scores to the `Pose.cache` dictionary: {tuple(_filtered_scores.keys())}. "
            + "Please use the `packed_pose = packed_pose.update_scores(**scores_dict)` syntax instead of "
            + "the `packed_pose.scores[key] = value` syntax to suppress this warning.",
            SyntaxWarning,
            stacklevel=2,
        )

    return packed_pose.update_scores(_filtered_scores)


@attr.define(kw_only=False, slots=True, frozen=True, auto_attribs=True)
class MessagePacking:
    """`MessagePack` base class for `PyRosettaCluster`."""

    pack: partial = attr.field(
        default=partial(msgpack.packb, use_bin_type=True),
        init=False,
        validator=attr.validators.instance_of(partial),
    )
    unpack: partial = attr.field(
        default=partial(msgpack.unpackb, raw=False),
        init=False,
        validator=attr.validators.instance_of(partial),
    )


@attr.define(kw_only=True, slots=False, frozen=False, auto_attribs=True)
class NonceCache:
    """Nonce cache base class for `PyRosettaCluster`."""

    instance_id: str = attr.field(
        validator=attr.validators.instance_of(str),
    )
    prk: Optional[bytes] = attr.field(
        default=None,
        repr=False,
        converter=lambda k: MaskedBytes(k) if isinstance(k, bytes) else k,
        validator=attr.validators.optional(attr.validators.instance_of(MaskedBytes)),
    )
    max_nonce: int = attr.field(
        validator=attr.validators.instance_of(int),
    )
    _seen: AbstractSet = attr.field(
        default=attr.Factory(set, takes_self=False),
        init=False,
        validator=attr.validators.instance_of(set),
    )
    _order: Deque[bytes] = attr.field(
        default=attr.Factory(
            lambda self: deque(maxlen=self.max_nonce),
            takes_self=True,
        ),
        init=False,
        validator=attr.validators.instance_of(deque),
    )
    _debug: bool = attr.field(
        default=False,
        init=False,
        validator=attr.validators.instance_of(bool),
    )

    def __attrs_pre_init__(self) -> None:
        """Pre-initialization hook for `NonceCache`."""
        MessagePacking.__init__(self)

    @staticmethod
    def _get_state(self) -> Dict[str, Any]:
        """
        A method used to override the default `NonceCache.__getstate__()` method that sets the value of the
        the 'prk' key (i.e., the pseudorandom key (PRK)) to `None` in the returned state.
        """
        state = self.__dict__.copy()
        state["prk"] = None

        return state

    @staticmethod
    def _on_worker() -> bool:
        """Test if the nonce cache is on a Dask worker."""
        try:
            _worker = get_worker()
            return True
        except Exception:
            return False

    def _cache_nonce(self, sealed: bytes) -> None:
        """
        Run Hash-based Message Authentication Code (HMAC) verification and cache nonces for replay protection
        without data decompression.
        """

        package = self.unpack(sealed)
        if not isinstance(package, dict) or package.get("v", None) != 1:
            _location = "Dask worker" if NonceCache._on_worker() else "head node"
            _err_msg = (
                f"Invalid sealed package or version on {_location} nonce cache! "
                + f"Received: {type(package)!r}"
            )
            raise ValueError(_err_msg)

        _instance_id = package["a"] # `str`: PyRosettaCluster instance identifier/App
        _data = package["d"] # `bytes`: Data bytestring
        _mac = package["m"] # `bytes`: MAC (HMAC tag)
        _nonce = package["n"] # `bytes`: Nonce
        _version = package["v"] # `int`: Version
        if _instance_id != self.instance_id:
            raise ValueError("`PyRosettaCluster` instance identifier mismatch in sealed package.")

        msg = self.pack([_instance_id, _data, _nonce, _version])
        _expected_mac = hmac_digest(bytes(self.prk), msg)
        if not compare_digest(_expected_mac, _mac):
            _location = "Dask worker" if NonceCache._on_worker() else "head node"
            _err_msg = (
                f"Task HMAC verification failed during nonce cache on {_location}!\n"
                + f"Expected: {_expected_mac!r}\n"
                + f"Value:    {_mac!r}\n"
            )
            raise SystemExit(_err_msg)

        if _nonce is not None:
            if _nonce in self._seen:
                # Replay protection
                _location = "Dask worker" if NonceCache._on_worker() else "head node"
                _err_msg = (
                    f"`PyRosettaCluster` detected a repeat nonce on the {_location} for the instance identifier "
                    + f"'{self.instance_id}', which might indicate a replay attack is in progress! "
                    + "Exiting process for security. Please ensure that `PyRosettaCluster(security=True)` "
                    + f"is enabled in future `PyRosettaCluster` simulations. Received: '{_nonce}'."
                )
                raise SystemExit(_err_msg)

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
                if NonceCache._on_worker():
                    print(f"Remote Dask worker ({get_worker().contact_address}) nonce cache: {_msg}")
                else:
                    print(f"Local head node nonce cache: {_msg}")


@attr.define(kw_only=True, slots=False, frozen=False, auto_attribs=True)
class Serialization:
    """Serialization base class for `PyRosettaCluster`."""

    instance_id: Optional[str] = attr.field(
        default=None,
        validator=attr.validators.optional(attr.validators.instance_of(str)),
    )
    prk: Optional[bytes] = attr.field(
        default=None,
        repr=False,
        converter=lambda k: MaskedBytes(k) if isinstance(k, bytes) else k,
        validator=attr.validators.optional(attr.validators.instance_of(MaskedBytes)),
    )
    compression: Optional[Union[str, bool]] = attr.field(
        default="xz",
        validator=attr.validators.optional(attr.validators.instance_of((str, bool))),
        converter=_parse_compression,
    )
    with_nonce: bool = attr.field(
        default=False,
        validator=attr.validators.instance_of(bool),
    )
    _nonce_size: int = attr.field(
        default=32, # bytes
        init=False,
        validator=attr.validators.instance_of(int),
    )

    def __attrs_pre_init__(self) -> None:
        """Pre-initialization hook for `Serialization`."""
        MessagePacking.__init__(self)

    def __attrs_post_init__(self) -> None:
        """Post-initialization hook for `Serialization`."""

        if self.compression == "xz":
            if "lzma" not in sys.modules:
                raise ImportError(
                    (
                        "Using 'xz' for compression requires installing the 'xz' package into "
                        + "your virtual environment. For installation instructions, visit:\n"
                        + "https://anaconda.org/anaconda/xz\n"
                        + "https://pypi.org/project/python-xz\n"
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
        A method used to override the default `Serialization.__getstate__()` method that sets  the value of the
        the 'prk' key (i.e., the pseudorandom key (PRK)) to `None` in the returned state.
        """
        state = self.__dict__.copy()
        state["prk"] = None

        return state

    @classmethod
    def zlib_compress(cls, obj: bytes) -> bytes:
        """Compress an object with `zlib` level 9."""
        return zlib.compress(obj, 9)

    def with_update_scores(func: CallableType) -> CallableType:
        """
        Decorator that caches detached `PackedPose.scores` items into the `Pose.cache` dictionary.

        *Warning*: This method uses the `pickle` module to deserialize pickled `Pose` objects. Using the
        `pickle` module is not secure, so please only run with input files you trust. Learn more about the
        `pickle` module and its security `here <https://docs.python.org/3/library/pickle.html>`_.
        """

        @wraps(func)
        def wrapper(self, obj: Any) -> Any:
            """Wrapper function to the `with_update_scores` decorator."""

            if isinstance(obj, PackedPose):
                obj = update_scores(obj)
            return func(self, obj)

        return cast(CallableType, wrapper)

    def requires_compression(func: CallableType) -> CallableType:
        """
        Decorator testing if compression is enabled, and skips compression if it is disabled.
        """

        @wraps(func)
        def wrapper(self, obj: Any) -> Any:
            """Wrapper function to the `requires_compression` decorator."""

            if all(x is not None for x in (self.encoder, self.decoder)):
                return func(self, obj)
            else:
                logging.debug("Compression/decompression is disabled.")
                return obj

        return cast(CallableType, wrapper)

    def _seal(self, data: bytes) -> bytes:
        """Seal data with `MessagePack`."""

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
        Unseal data with `MessagePack` and perform Hash-based Message Authentication Code (HMAC) verification.
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

    @with_update_scores
    @requires_compression
    def compress_packed_pose(self, packed_pose: Optional[PackedPose]) -> Optional[bytes]:
        """
        Compress a `PackedPose` object with the custom serialization module. If the 'packed_pose' argument
        value is `None`, then just return `None`.

        Args:
            `packed_pose`:
                An input `PackedPose` object to compress. If `None`, then just return `None`.

        Returns:
            A `bytes` object representing the compressed `PackedPose` object, or `None`.

        Raises:
            `TypeError` if the `packed_pose` argument value is not of type `NoneType` or `PackedPose`.
        """

        if packed_pose is None:
            compressed_packed_pose = None
        elif isinstance(packed_pose, PackedPose):
            compressed_packed_pose = self.encoder(io.to_pickle(packed_pose))
        else:
            raise TypeError(
                "The `packed_pose` argument value must be of type `NoneType` or `PackedPose`. "
                + f"Received: {type(packed_pose)}"
            )

        return compressed_packed_pose

    @requires_compression
    def decompress_packed_pose(self, compressed_packed_pose: Optional[bytes]) -> Optional[PackedPose]:
        """
        Decompress a `bytes` object with the custom serialization module and secure implementation of the
        `pickle` module. If the `compressed_packed_pose` argument value is `None`, then just return `None`.

        *Warning*: This function uses the `pickle` module to deserialize pickled `Pose` objects. Using the
        `pickle` module is not secure, so please only run with input files you trust. Learn more about the
        `pickle` module and its security `here <https://docs.python.org/3/library/pickle.html>`_.

        Args:
            `compressed_packed_pose`:
                An input `bytes` object to decompress. If `None`, then just return `None`.

        Returns:
            A `PackedPose` object representing the decompressed `bytes` object, or `None`.

        Raises:
            `TypeError` if the `compressed_packed_pose` argument value is not of type `NoneType` or `bytes`.
        """

        if compressed_packed_pose is None:
            packed_pose = None
        elif isinstance(compressed_packed_pose, bytes):
            packed_pose = io.to_packed(SecureSerializerBase.secure_loads(self.decoder(compressed_packed_pose)))
        else:
            raise TypeError(
                "The `compressed_packed_pose` argument value must be of type `NoneType` or `bytes`. "
                + f"Received: {type(compressed_packed_pose)}"
            )

        return packed_pose

    def loads_object(self, compressed_obj: bytes) -> Any:
        """
        Unseal data and run the `cloudpickle.loads` method.

        *Warning*: This method uses the `cloudpickle` module to deserialize objects. Using the `cloudpickle`
        module is not secure, so please only run with input data you trust. Learn more about the `cloudpickle`
        module and its security `here <https://github.com/cloudpipe/cloudpickle>`_ and
        `here <https://docs.python.org/3/library/pickle.html>`_.
        """

        data = self._unseal(compressed_obj)
        buffer = self.decoder(data) if self.decoder else data

        return cloudpickle.loads(buffer)

    def dumps_object(self, obj: Any) -> bytes:
        """Run the `cloudpickle.dumps` method and seal data."""

        pickled = cloudpickle.dumps(obj)
        buffer = self.encoder(pickled) if self.encoder else pickled

        return self._seal(buffer)

    def compress_kwargs(self, kwargs: Dict[str, Any]) -> bytes:
        """
        Compress a `dict` object with the `cloudpickle` and custom serialization modules.

        Args:
            `kwargs`:
                An input `dict` object to compress.

        Returns:
            A `bytes` object representing the compressed `dict` object.

        Raises:
            `TypeError` if the 'kwargs' argument value is not of type `dict`.
        """

        if isinstance(kwargs, dict):
            return self.dumps_object(kwargs)
        else:
            raise TypeError(
                "The `kwargs` argument value must be of type `dict`. "
                + f"Received: {type(kwargs)}"
            )

    def decompress_kwargs(self, compressed_kwargs: bytes) -> Dict[str, Any]:
        """
        Decompress a `bytes` object with the custom serialization and `cloudpickle` modules.

        *Warning*: This method uses the `cloudpickle` module to deserialize objects. Using the `cloudpickle`
        module is not secure, so please only run with input data you trust. Learn more about the `cloudpickle`
        module and its security `here <https://github.com/cloudpipe/cloudpickle>`_ and
        `here <https://docs.python.org/3/library/pickle.html>`_.

        Args:
            `compressed_kwargs`:
                An input `bytes` object to decompress.

        Returns:
            A `dict` object representing the decompressed `bytes` object.

        Raises:
            `TypeError` if the `compressed_packed_pose` argument value is not of type `bytes`.
            `TypeError` if the decompressed object is not of type `dict`.
        """

        if isinstance(compressed_kwargs, bytes):
            kwargs = self.loads_object(compressed_kwargs)
            if not isinstance(kwargs, dict):
                raise TypeError(f"Decompressed object must be of type `dict`. Received: {type(kwargs)}")
            return kwargs
        else:
            raise TypeError(
                "The `compressed_kwargs` argument value must be of type `bytes`. "
                + f"Received: {type(compressed_kwargs)}"
            )

    def compress_object(self, obj: Any) -> bytes:
        """
        Compress an object with the `cloudpickle` and custom serialization modules.

        Args:
            `obj`:
                An input object to compress.

        Returns:
            A `bytes` object representing the compressed object.
        """
        return self.dumps_object(obj)

    def decompress_object(self, compressed_obj: bytes) -> Any:
        """
        Decompress a `bytes` object with the custom serialization and `cloudpickle` modules.

        *Warning*: This method uses the `cloudpickle` module to deserialize objects. Using the `cloudpickle`
        module is not secure, so please only run with input data you trust. Learn more about the `cloudpickle`
        module and its security `here <https://github.com/cloudpipe/cloudpickle>`_ and
        `here <https://docs.python.org/3/library/pickle.html>`_.

        Args:
            `compressed_obj`:
                An input `bytes` object to decompress.

        Returns:
            An object representing the decompressed `bytes` object.

        Raises:
            `TypeError` if the `compressed_obj` argument value is not of type `bytes`.
        """

        if isinstance(compressed_obj, bytes):
            return self.loads_object(compressed_obj)
        else:
            raise TypeError(
                "The `compressed_obj` argument value must be of type `bytes`. "
                + f"Received: {type(compressed_obj)}"
            )

    @classmethod
    def deepcopy_kwargs(cls, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """
        The `cloudpickle` module makes it possible to serialize Python constructs not supported by the default
        `pickle` module from the Python standard library.

        *Warning*: This method uses the `cloudpickle` module to serialize and subsequently deserialize `dict`
        objects. Using the `cloudpickle` module is not secure, so please only run with input data you trust.
        Learn more about the `cloudpickle` module and its security
        `here <https://github.com/cloudpipe/cloudpickle>`_ and
        `here <https://docs.python.org/3/library/pickle.html>`_.

        Args:
            `kwargs`:
                A `dict` object to be deep copied.

        Returns:
            A deep copy of the `dict` object.

        Raises:
            `TypeError` if the `kwargs` argument value is not of type `dict`.
        """

        if isinstance(kwargs, dict):
            return cloudpickle.loads(cloudpickle.dumps(kwargs))
        else:
            raise TypeError(
                "The `kwargs` argument value must be of type `dict`. "
                + f"Received: {type(kwargs)}"
            )
