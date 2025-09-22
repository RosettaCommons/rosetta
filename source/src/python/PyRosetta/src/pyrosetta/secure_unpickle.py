# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.
# Secure unpickling in PyRosetta


__author__ = "Jason C. Klima"


import base64
import hashlib
import hmac
import importlib
import io
import pickle
import pyrosetta.rosetta  # noqa
import sys

from pathlib import Path
from functools import lru_cache, partial
from typing import AbstractSet, Any, Dict, Iterable, NoReturn, Optional, Tuple, Union


# Default secure python builtins:
SECURE_PYTHON_BUILTINS: AbstractSet[str] = {
    "bool",
    "bytearray",
    "bytes",
    "complex",
    "dict",
    "float",
    "frozenset",
    "int",
    "list",
    "NoneType",
    "range",
    "set",
    "slice",
    "str",
    "tuple",
    "type",
}

# Default secure packages, not including PyRosetta:
SECURE_EXTRA_PACKAGES: Tuple[str, ...] = ("numpy",)

# Modules that will not be callable targets via the `pickle` module:
BLOCKED_PACKAGES: AbstractSet[str] = {
    "subprocess",  # Block `subprocess.*` calls
    "ctypes",  # Block arbitrary code
    # Block pickle and relatives:
    "cloudpickle",
    "dill",
    "joblib",
    "pickle",
    "zodbpickle",
    # Block process spawning:
    "billiard",
    "celery",
    "dask",
    "distributed",
    "multiprocessing",
    # Other security concerns <https://docs.python.org/3/library/security_warnings.html>:
    "base64",
    "hashlib",
    "http",
    "logging",
    "random",
    "shelve",
    "ssl",
    "tempfile",
    "xml",
    "zipfile",
}

# Globals that will always be blocked:
BLOCKED_GLOBALS: AbstractSet[Tuple[str, str]] = {
    ("builtins",  "__import__"),
    ("builtins",  "compile"),
    ("builtins",  "eval"),
    ("builtins",  "exec"),
    ("builtins",  "open"),
    ("http",      "server"),
    ("importlib", "import_module"),
    ("os",        "_exit"),
    ("os",        "popen"),
    ("os",        "system"),
    ("posix",     "popen"),
    ("posix",     "system"),
    ("sys",       "exit"),
}

# Package:prefix pairs for specific method prefixes that will always be blocked:
BLOCKED_PREFIXES: Dict[str, Tuple[str, ...]] = {
    "os": (
        "execl",        "execle",       "execlp",      "execlpe",
        "execv",        "execve",       "execvp",      "execvpe",
        "spawnl",       "spawnle",      "spawnlp",     "spawnlpe",
        "spawnv",       "spawnve",      "spawnvp",     "spawnvpe",
        "abort",        "access",       "chdir",       "chflags",
        "chmod",        "chown",        "chroot",      "close",
        "fchdir",       "fchmod",       "fchown",      "fdopen",
        "fork",         "forkpty",      "fsync",       "kill",
        "lchmod",       "lchown",       "link",        "pipe",
        "posix_spawn",  "posix_spawnp", "remove",      "rename",
        "rmdir",        "startfile",    "sync",        "sys",
        "unlink",       "write",
    ),
    "posix": (
        "close",        "execv",        "execve",      "kill",
        "open",
    ),
    "shutil": (
        "chown",        "make_archive", "move",        "rmtree",
    ),
    "sys": (
        "addaudithook", "setprofile",   "settrace",
    ),
}


# Hash-based Message Authentication Code (HMAC) key for data integrity:
HASHMOD: partial = partial(hashlib.blake2s, digest_size=16, salt=b'cache')
HMAC_SIZE: int = HASHMOD().digest_size
HMAC_KEY: Optional[bytes] = None  # Disabled by default

def set_unpickle_hmac_key(key: Optional[bytes]) -> None:
    """
    Set the global Hash-based Message Authentication Code (HMAC) key
    for `Pose.cache` score object secure serialization.
    """
    global HMAC_KEY
    if key is not None and not isinstance(key, bytes):
        raise TypeError(
            "The 'key' argument parameter must be a `bytes` or `NoneType` object. "
            + "Received: %s" % type(key)
        )
    HMAC_KEY = key

def get_unpickle_hmac_key() -> Optional[bytes]:
    """
    Get the global Hash-based Message Authentication Code (HMAC) key
    for `Pose.cache` score object secure serialization.
    """
    return HMAC_KEY


# `UnpicklingError` exception subclasses:

class UnpickleCompatibilityError(pickle.UnpicklingError):
    """
    Subclass of `pickle.UnpicklingError` raised when an unpickle-allowed module
    cannot be resolved due to a Python package version or environment mismatch
    from that used to pickle the module.
    """
    def __init__(self, module: str, name: str) -> None:
        _top_package = _split_top_package(module)
        _msg = (
            "Unable to unpickle the allowed '%s.%s' symbol due to " % (module, name,)
            + "a Python version mismatch, a virtual environment mismatch, or a "
            + "missing dependency of the pickled data. Please install or upgrade the "
            + "required package and try again: %r" % (_top_package,)
        )
        super().__init__(_msg)
        self.module = module
        self.name = name
        self._top_package = _top_package

class UnpickleIntegrityError(pickle.UnpicklingError):
    """Subclass of `pickle.UnpicklingError` raised on failed HMAC verification."""
    def __init__(self, *args: Any) -> None:
        super().__init__(*args)

class UnpickleSecurityError(pickle.UnpicklingError):
    """
    Subclass of `pickle.UnpicklingError` raised when pickled objects
    reference disallowed globals and modules.
    """
    def __init__(self, module: str, name: str, allowed: Tuple[str, ...]) -> None:
        _top_package = _split_top_package(module)
        _allowed = tuple(sorted(set(("pyrosetta",) + allowed)))
        _disallowed = get_disallowed_packages()
        _msg = (
            "Disallowed unpickling of the '%s.%s' namespace!\n" % (module, name,)
            + "The currently allowed packages to be securely unpickled are: %s\n" % (_allowed,)
            + "The received object requires the %r package. " % (_top_package,)
        )
        if (
            module in _disallowed
            or "%s.%s" % (module, name,) in _disallowed
            or "%s.%s*" % (module, name,) in _disallowed
        ):
            _msg += (
                "However, the '%s.%s' namespace cannot be added to the set of trusted packages " % (module, name,)
                + "since it is permanently disallowed! To view the set of permanently disallowed packages, "
                + "please run:\n    `pyrosetta.get_disallowed_packages()`\n"
            )
        elif (module == "builtins" and name not in SECURE_PYTHON_BUILTINS):
            _msg += (
                "However, the '%s.%s' method cannot be added to the set of trusted packages " % (module, name,)
                + "since it is not allowed! Please consider reporting an issue if the %r python builtins " % (name,)
                + "needs to be unpickled for your application."
            )
        else:
            if _top_package in allowed:
                _msg += (
                    "However, the %r package is already a trusted package, so the %r module could not be resolved! " % (_top_package, module,)
                    + "Please consider reporting an issue if the %r module needs to be unpickled for your application." % (module,)
                )
            else:
                _msg += (
                    "To add it to the set of trusted packages, please run the following then try again:\n"
                    + "    `pyrosetta.add_secure_package(%r)`\n" % (_top_package,)
                )
        super().__init__(_msg)
        self.module = module
        self.name = name
        self.allowed = allowed
        self._top_package = _top_package
        self._allowed = _allowed
        self._disallowed = _disallowed


# Methods to update the unpickle-allowed list of secure packages:

def add_secure_package(package: str) -> None:
    """
    Add a secure package by top-level name to the unpickle-allowed list.
    """
    if not package:
        return None
    _top_package = _split_top_package(str(package))
    _secure_packages = list(get_secure_packages())
    if _top_package not in _secure_packages:
        set_secure_packages(tuple(_secure_packages + [_top_package]))

def clear_secure_packages() -> None:
    """
    Remove all secure packages, excluding 'pyrosetta' which is always implicitly allowed.
    """
    set_secure_packages(tuple())

def get_secure_packages() -> Tuple[str, ...]:
    """
    Return the extra secure packages currently allowed, excluding 'pyrosetta' which is
    always implicitly allowed.
    """
    return SECURE_EXTRA_PACKAGES

def remove_secure_package(package: str) -> None:
    """
    Remove a secure package by top-level name if present in the unpickle-allowed list.
    """
    if not package:
        return None
    _top_package = _split_top_package(str(package))
    _secure_packages = [p for p in get_secure_packages() if p != _top_package]
    set_secure_packages(tuple(_secure_packages))

def set_secure_packages(packages: Iterable[str]) -> None:
    """
    Set the secure extra packages in the unpickle-allowed list, excluding 'pyrosetta' which is
    always implicitly allowed.

    Example:
        `set_secure_packages(('numpy', 'pandas'))`
    """
    global SECURE_EXTRA_PACKAGES

    if not isinstance(packages, (list, tuple, set)):
        raise TypeError(
            "The 'packages' argument parameter must be a `list`, `tuple`, or `set` object. "
            + "Received: %s" % type(packages)
        )
    _seen = set()
    _out = []
    for package in packages:
        if not isinstance(package, str):
            raise TypeError(
                f"The 'packages' argument parameter items must be of type `str`. Received: {type(package)}"
            )
        if not package:
            continue
        _top_package = _split_top_package(package)
        if _top_package not in _seen:
            _seen.add(_top_package)
            _out.append(_top_package)
    SECURE_EXTRA_PACKAGES = tuple(sorted(_out))

def get_disallowed_packages() -> Tuple[str, ...]:
    """
    Return a `tuple` of packages and methods that are permanently disallowed
    from being unpickled in PyRosetta, where '*' matches any string.
    """
    disallowed = set()
    for _package in BLOCKED_PACKAGES:
        disallowed.add(_package)
    for _package, _method in BLOCKED_GLOBALS:
        disallowed.add("%s.%s" % (_package, _method,))
    for _package, _methods in BLOCKED_PREFIXES.items():
        for _method in _methods:
            disallowed.add("%s.%s*" % (_package, _method,))

    return tuple(sorted(disallowed))

def _split_top_package(module: str) -> str:
    return module.split(".", 1)[0]


# Unpickle methods:

class ModuleCache(object):
    """
    Resolve modules and packages by path, and determine if they are allowed or disallowed.
    """
    @staticmethod
    @lru_cache(maxsize=1)
    def _rosetta_module() -> object:
        _module = sys.modules.get("pyrosetta.rosetta", None)
        if _module is None:
            __import__("pyrosetta.rosetta")
            _module = sys.modules.get("pyrosetta.rosetta", None)
            if _module is None:
                raise ImportError("pyrosetta.rosetta")

        return _module

    @staticmethod
    @lru_cache(maxsize=1)
    def _rosetta_origin() -> Optional[Path]:
        _rosetta_module = ModuleCache._rosetta_module()
        _rosetta_spec = getattr(_rosetta_module, "__spec__", None)
        _rosetta_origin = getattr(_rosetta_spec, "origin", None) or getattr(_rosetta_module, "__file__", None)
        if _rosetta_origin:
            _rosetta_origin = Path(_rosetta_origin).resolve()

        return _rosetta_origin

    @staticmethod
    @lru_cache(maxsize=1024, typed=True)
    def _package_base_dir(package_name: str) -> Optional[Path]:
        try:
            _package = importlib.import_module(package_name)
        except Exception:
            return None
        _file = getattr(_package, "__file__", None)

        return Path(_file).resolve().parent if _file else None

    @staticmethod
    @lru_cache(maxsize=4096, typed=True)
    def _module_file(module_name: str) -> Optional[Path]:
        try:
            _module = importlib.import_module(module_name)
        except Exception:
            return None
        _file = getattr(_module, "__file__", None)

        return Path(_file).resolve() if _file else None

    @staticmethod
    def _is_relative_to(path: Path, base: Path) -> bool:
        try:
            path.relative_to(base)
            return True
        except Exception:
            return False

    @staticmethod
    def _is_under_package(module: str, package: str) -> bool:
        _base_dir = ModuleCache._package_base_dir(package)
        _module_file = ModuleCache._module_file(module)

        return bool(
            _base_dir
            and _module_file
            and ModuleCache._is_relative_to(_module_file, _base_dir)
        )

    @staticmethod
    def _is_under_rosetta(module: str) -> bool:
        if not (module == "pyrosetta.rosetta" or module.startswith("pyrosetta.rosetta.")):
            return False
        # Check if submodule has an origin identical to the 'pyrosetta.rosetta' origin
        _module = sys.modules.get(module, None)
        if _module is not None:
            _module_spec = getattr(_module, "__spec__", None)
            _module_origin = getattr(_module_spec, "origin", None) or getattr(_module, "__file__", None)
            if _module_origin:
                return Path(_module_origin).resolve() == ModuleCache._rosetta_origin()
        # Otherwise, walk down attributes of imported virtual submodule
        if ModuleCache._walk_rosetta_module(module) is None:
            return False
        # Otherwise, attribute exists under the 'pyrosetta.rosetta' module
        return True

    @staticmethod
    def _walk_rosetta_module(module: str) -> Any:
        assert (module == "pyrosetta.rosetta" or module.startswith("pyrosetta.rosetta."))
        _obj = ModuleCache._rosetta_module()
        for _name in module.split(".")[2:]:  # Skip 'pyrosetta.rosetta'
            _obj = getattr(_obj, _name, None)
            if _obj is None:
                return None

        return _obj

    @staticmethod
    def _is_allowed_module(module: str) -> bool:
        # Always trust PyRosetta modules by path
        if (module == "pyrosetta" or module.startswith("pyrosetta.")):
            if (module == "pyrosetta.rosetta" or module.startswith("pyrosetta.rosetta.")):
                return ModuleCache._is_under_rosetta(module)
            return ModuleCache._is_under_package(module, "pyrosetta")
        else: # Maybe trust other modules
            _top_package = _split_top_package(module)
            if _top_package in SECURE_EXTRA_PACKAGES and ModuleCache._is_under_package(module, _top_package):
                return True
            else:
                return False

    @staticmethod
    def _get_allowed_module_attr(module: str, name: str) -> Any:
        if (module == "pyrosetta.rosetta" or module.startswith("pyrosetta.rosetta.")):
            # Prevent re-import; instead walk down attributes of imported virtual submodule
            _module = ModuleCache._walk_rosetta_module(module)
            if _module is None:
                raise UnpickleCompatibilityError(module, name)
        elif module == "builtins":
            if name == "NoneType":
                return type(None)
            _module = sys.modules["builtins"]
        else:
            # Prevent re-import if the module is already imported
            try:
                _module = sys.modules.get(module, None) or importlib.import_module(module)
            except ImportError as ex:
                raise UnpickleCompatibilityError(module, name) from ex
        try:
            return getattr(_module, name)
        except AttributeError as ex:
            raise UnpickleCompatibilityError(module, name) from ex


class SecureUnpickler(pickle.Unpickler):
    """
    Secure subclass of `pickle.Unpickler` predicated on allowed and disallowed globals, modules, and prefixes.
    """
    def __init__(self, file: io.BytesIO, *, stream_protocol: int = -1) -> None:
        super().__init__(file)
        if not isinstance(stream_protocol, int):
            raise TypeError(
                "The 'stream_protocol' keyword argument parameter must be "
                + "an `int` object. Received: %s" % (type(stream_protocol),)
            )
        self._stream_protocol = stream_protocol

    def find_class(self, module: str, name: str) -> Union[Any, NoReturn]:
        if module in BLOCKED_PACKAGES:
            raise UnpickleSecurityError(module, name, get_secure_packages())
        if (module, name) in BLOCKED_GLOBALS:
            raise UnpickleSecurityError(module, name, get_secure_packages())
        if module in BLOCKED_PREFIXES and any(name.startswith(prefix) for prefix in BLOCKED_PREFIXES[module]):
            raise UnpickleSecurityError(module, name, get_secure_packages())
        # Builtins:
        if module == "builtins":
            if name in SECURE_PYTHON_BUILTINS:
                return ModuleCache._get_allowed_module_attr(module, name)
            raise UnpickleSecurityError(module, name, get_secure_packages())
        # Maybe include `copyreg` unpickle helper functions, depending on incoming stream protocol:
        if module == "copyreg":
            if (0 <= self._stream_protocol <= 1 and name == "_reconstructor") or (
                2 <= self._stream_protocol <= pickle.HIGHEST_PROTOCOL and name in ("__newobj__", "__newobj_ex__",)
            ):
                return ModuleCache._get_allowed_module_attr(module, name)
            raise UnpickleSecurityError(module, name, get_secure_packages())
        if ModuleCache._is_allowed_module(module):
            return ModuleCache._get_allowed_module_attr(module, name)

        raise UnpickleSecurityError(module, name, get_secure_packages())


# Secure serialization for `PackedPose` and `Pose.cache` score objects:

class SecureSerializerBase(object):
    """
    Base class for `PackedPose`, `Pose`, and `Pose.cache` score
    object secure serialization.
    """
    _encoder: str = "utf-8"
    _pickle_protocol: int = pickle.DEFAULT_PROTOCOL

    @staticmethod
    def to_pickle(value: Any) -> Union[bytes, NoReturn]:
        try:
            return pickle.dumps(value, protocol=SecureSerializerBase._pickle_protocol)
        except (TypeError, OverflowError, MemoryError, pickle.PicklingError) as ex:
            raise Exception(
                "%s: %s. Only pickle-serializable object types are " % (type(ex).__name__, ex,)
                + "allowed to be dumped. Received: %r. %s" % (type(value), ex,)
            ) from ex

    @staticmethod
    def from_base64(value: Union[str, bytes]) -> bytes:
        return base64.b64decode(value, validate=True)

    @staticmethod
    def to_base64(value: bytes) -> str:
        return base64.b64encode(value).decode(SecureSerializerBase._encoder)

    @staticmethod
    def secure_loads(value: bytes) -> Union[Any, NoReturn]:
        """Secure replacement for `pickle.loads`."""
        stream = io.BytesIO(value)
        stream_protocol = SecureSerializerBase._get_stream_protocol(value)
        try:
            return SecureUnpickler(stream, stream_protocol=stream_protocol).load()
        except (UnpickleSecurityError, UnpickleIntegrityError, UnpickleCompatibilityError):
            raise
        except MemoryError:
            raise
        except KeyboardInterrupt:
            raise
        except pickle.UnpicklingError as ex:
            raise pickle.UnpicklingError(
                f"{ex}. PyRosetta secure unpickle failed with stream protocol {stream_protocol}."
            ) from ex
        except (TypeError, OverflowError, EOFError) as ex:
            raise Exception(
                "%s: %s. Could not deserialize value of type %r." % (type(ex).__name__, ex, type(value),)
            ) from ex

    @staticmethod
    def secure_from_base64_pickle(string: str) -> Any:
        raw = SecureSerializerBase.from_base64(string)
        key = get_unpickle_hmac_key()
        if key is not None:
            raw = SecureSerializerBase._verify_and_remove_hmac_tag(key, raw)

        return SecureSerializerBase.secure_loads(raw)

    @staticmethod
    def secure_to_base64_pickle(obj: Any) -> str:
        obj = SecureSerializerBase.to_pickle(obj)
        key = get_unpickle_hmac_key()
        if key is not None:
            obj = SecureSerializerBase._prepend_hmac_tag(key, obj)

        return SecureSerializerBase.to_base64(obj)

    @staticmethod
    def _get_hmac_tag(key: bytes, data: bytes) -> bytes:
        return hmac.digest(key, data, HASHMOD)

    @staticmethod
    def _prepend_hmac_tag(key: bytes, data: bytes) -> bytes:
        return SecureSerializerBase._get_hmac_tag(key, data) + data

    @staticmethod
    def _verify_and_remove_hmac_tag(key: bytes, signed_data: bytes) -> Union[bytes, NoReturn]:
        if len(signed_data) < HMAC_SIZE:
            raise UnpickleIntegrityError("Signed object is too short.")
        hmac_tag, data = signed_data[:HMAC_SIZE], signed_data[HMAC_SIZE:]
        expected_hmac_tag = SecureSerializerBase._get_hmac_tag(key, data)
        if not hmac.compare_digest(hmac_tag, expected_hmac_tag):
            raise UnpickleIntegrityError(
                "PyRosetta secure unpickler HMAC verification failed! The HMAC key may have "
                + "changed or the data has been tampered with after the data was pickled."
            )

        return data

    @staticmethod
    def _get_stream_protocol(obj: bytes) -> int:
        if len(obj) >= 2 and obj[0] == pickle.PROTO[0]:
            return obj[1] if 0 <= obj[1] <= pickle.HIGHEST_PROTOCOL else -1
        else:
            return -1
