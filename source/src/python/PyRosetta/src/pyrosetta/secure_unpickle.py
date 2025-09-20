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
    "memoryview",
    "set",
    "str",
    "tuple",
}

# Default secure packages, not including PyRosetta:
SECURE_EXTRA_PACKAGES: Tuple[str, ...] = ("numpy",)

# Modules that will not be callable targets via the `pickle` module:
BLOCKED_PACKAGES: AbstractSet[str] = {
    "subprocess",  # Block `subprocess.*` calls (including `run`, `check_output`, `Popen`, etc.)
    "ctypes",  # Block arbitrary code
    # Block pickle and alternatives:
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
}

# Globals that will always be blocked:
BLOCKED_GLOBALS: AbstractSet[Tuple[str, str]] = {
    ("builtins", "eval"),
    ("builtins", "exec"),
    ("builtins", "compile"),
    ("builtins", "__import__"),
    ("builtins", "open"),
    ("os", "system"),
    ("os", "popen"),
    ("posix", "system"),
    ("posix", "popen"),
    ("sys", "exit"),
    ("importlib", "import_module")
}

# Package:prefix pairs for specific method prefixes that will always be blocked:
BLOCKED_PREFIXES: Dict[str, Tuple[str, ...]] = {
    "os": (
        "execl", "execle", "execlp", "execlpe",
        "execv", "execve", "execvp", "execvpe",
        "spawnl", "spawnle", "spawnlp", "spawnlpe",
        "spawnv", "spawnve", "spawnvp", "spawnvpe",
        "posix_spawn", "posix_spawnp",
        "remove", "rmdir", "unlink",
        "startfile",
    ),
    "posix": ("close", "execv", "execve", "kill", "open",),
    "shutil": ("chown", "make_archive", "move", "rmtree",),
    "sys": ("addaudithook", "setprofile", "settrace",),
}


# Hash-based Message Authentication Code (HMAC) key for data integrity:
HASHMOD: partial = partial(hashlib.blake2s, digest_size=8)
HMAC_SIZE: int = HASHMOD().digest_size
HMAC_KEY: Optional[bytes] = None  # Disabled by default

def set_unpickle_hmac_key(key: Optional[bytes]) -> None:
    """
    Set the global Hash-based Message Authentication Code (HMAC) key
    for Pose.cache` score object secure serialization.
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
    for Pose.cache` score object secure serialization.
    """
    return HMAC_KEY


# `UnpicklingError` exception subclasses:

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
        self.module = module
        self.name = name
        self.allowed = allowed
        top_package = _split_top_package(module)
        allowed_sorted = tuple(sorted(set(("pyrosetta",) + self.allowed)))
        msg = (
            "Disallowed unpickling of the '%s.%s' namespace!\n" % (module, name)
            + "The currently allowed packages to be securely unpickled are: %s\n" % (allowed_sorted,)
            + "The received object requires the %r package. " % (top_package,)
        )
        disallowed = get_disallowed_packages()
        if (
            module in disallowed
            or "%s.%s" % (module, name,) in disallowed
            or "%s.%s*" % (module, name,) in disallowed
        ):
            msg += (
                "However, the '%s.%s' namespace cannot be added to the set of trusted packages " % (module, name)
                + "since it is permanently disallowed! To view the set of permanently disallowed packages, "
                + "please run:\n    `pyrosetta.get_disallowed_packages()`\n"
            )
        else:
            msg += (
                "To add it to the set of trusted packages, please run the following then try again:\n"
                + "    `pyrosetta.add_secure_package(%r)`\n" % (top_package,)
            )
        super().__init__(msg)


# Methods to update the unpickle-allowed list of secure packages:

def add_secure_package(package: str) -> None:
    """
    Add a secure package by top-level name to the unpickle-allowed list.
    """
    if not package:
        return None
    top = _split_top_package(str(package))
    current = list(get_secure_packages())
    if top not in current:
        set_secure_packages(tuple(current + [top]))

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
    top = _split_top_package(str(package))
    current = [p for p in get_secure_packages() if p != top]
    set_secure_packages(tuple(current))

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
    seen = set()
    out = []
    for package in packages:
        if not package:
            continue
        top = _split_top_package(str(package))
        if top not in seen:
            seen.add(top)
            out.append(top)
    SECURE_EXTRA_PACKAGES = tuple(out)

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
    Resolve modules and packages by path, and determine if they are allowed or blocked.
    """
    @staticmethod
    @lru_cache(maxsize=None)
    def _package_base_dir(package_name: str) -> Optional[Path]:
        try:
            package = importlib.import_module(package_name)
        except Exception:
            return None
        f = getattr(package, "__file__", None)
        return Path(f).resolve().parent if f else None

    @staticmethod
    @lru_cache(maxsize=None)
    def _module_file(module: str) -> Optional[Path]:
        try:
            mod = importlib.import_module(module)
        except Exception:
            return None
        f = getattr(mod, "__file__", None)
        return Path(f).resolve() if f else None

    @staticmethod
    def _is_relative_to(p: Path, base: Path) -> bool:
        try:
            p.relative_to(base)
            return True
        except Exception:
            return False

    @staticmethod
    def _is_under_package(module: str, package: str) -> bool:
        base = ModuleCache._package_base_dir(package)
        mf = ModuleCache._module_file(module)
        return bool(base and mf and ModuleCache._is_relative_to(mf, base))

    @staticmethod
    def _is_allowed_module(module: str) -> bool:
        # Always trust PyRosetta modules by path
        if (module == "pyrosetta"):
            return True
        if module.startswith("pyrosetta."):
            if (
                (module == "pyrosetta.rosetta")
                or module.startswith("pyrosetta.rosetta.")
                or ModuleCache._is_under_package(module, "pyrosetta")
            ):
                return True
            else:
                raise ModuleNotFoundError(module)
        else: # Maybe trust other modules
            top = _split_top_package(module)
            if top in SECURE_EXTRA_PACKAGES and ModuleCache._is_under_package(module, top):
                return True
            else:
                return False


class SecureUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module in BLOCKED_PACKAGES:
            raise UnpickleSecurityError(module, name, get_secure_packages())
        if (module, name) in BLOCKED_GLOBALS:
            raise UnpickleSecurityError(module, name, get_secure_packages())
        if module in BLOCKED_PREFIXES and any(name.startswith(prefix) for prefix in BLOCKED_PREFIXES[module]):
            raise UnpickleSecurityError(module, name, get_secure_packages())
        # Builtins:
        if module == "builtins" and name in SECURE_PYTHON_BUILTINS:
            return getattr(sys.modules["builtins"], name)
        # For pickle protocols 0 and 1, include `copyreg` unpickle helper function:
        if SecureSerializerBase._pickle_protocol in (0, 1):
            if module == "copyreg" and name == "_reconstructor":
                __import__(module)
                return getattr(sys.modules[module], name)
        if ModuleCache._is_allowed_module(module):
            __import__(module)
            return getattr(sys.modules[module], name) 

        raise UnpickleSecurityError(module, name, get_secure_packages())


# Secure serialization for `PackedPose` and `Pose.cache` score objects:

class SecureSerializerBase(object):
    """
    Base class for for `PackedPose`, `Pose`, and `Pose.cache` score
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
            )

    @staticmethod
    def from_base64(value: str) -> bytes:
        return base64.b64decode(value, validate=True)

    @staticmethod
    def to_base64(value: bytes) -> str:
        return base64.b64encode(value).decode(SecureSerializerBase._encoder)

    @staticmethod
    def secure_loads(value: bytes) -> Any:
        """Secure replacement for `pickle.loads`."""
        try:
            return SecureUnpickler(io.BytesIO(value)).load()
        except (TypeError, OverflowError, MemoryError, EOFError) as ex:
            raise Exception(
                "%s: %s. Could not deserialize value of type %r." % (type(ex).__name__, ex, type(value),)
            )
        except UnpickleSecurityError as ex:
            raise UnpickleSecurityError(ex.module, ex.name, ex.allowed) from ex

    @staticmethod
    def secure_from_base64_pickle(string: str) -> Any:
        raw = SecureSerializerBase.from_base64(
            string.encode(SecureSerializerBase._encoder)
        )
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
