# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


try:
    import msgpack
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.hkdf' requires the "
        + "third-party package 'msgpack' as a dependency!\n"
        + "Please install the package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/msgpack/\n"
    )
    raise

import hashlib
import hmac
import types

from typing import Optional


HASHMOD: types.BuiltinFunctionType = hashlib.sha256
DERIVED_KEY_LEN: int = 32 


def compare_digest(a: bytes, b: bytes) -> bool:
    """Return `a == b` result."""
    return hmac.compare_digest(a, b)


def hmac_digest(key: bytes, data: bytes) -> bytes:
    """Return the digest of a new hash-based message authentication code (HMAC) object."""
    return hmac.new(key, data, HASHMOD).digest()


def hkdf_extract(ikm: bytes, salt: Optional[bytes]) -> bytes:
    """
    Extract method for hash-based message authentication code (HMAC)-based key derivation
    function (HKDF), taking input keying material (IKM) and optional salt.
    """
    if salt is None:
        salt = b"\x00" * HASHMOD().digest_size

    return hmac_digest(salt, ikm)


def hkdf_expand(prk: bytes, info: bytes) -> bytes:
    """
    Expand method for hash-based message authentication code (HMAC)-based key derivation
    function (HKDF), taking a pseudorandom key (PRK) and application-specific info.
    """
    okm = b"" 
    prev_block = b""
    hash_len = HASHMOD().digest_size
    blocks = (DERIVED_KEY_LEN + hash_len - 1) // hash_len
    for counter in range(1, blocks + 1):
        data = prev_block + info + bytes([counter])
        block = hmac_digest(prk, data)
        okm += block
        prev_block = block

    return okm[:DERIVED_KEY_LEN]


def derive_task_key(passkey: bytes, task_id: str) -> bytes:
    """
    Derive a per-task secret key using the extract-and-expand hash-based message
    authentication code (HMAC)-based key derivation function (HKDF).
    """
    info = msgpack.packb(["PyRosettaCluster", task_id], use_bin_type=True)
    prk = hkdf_extract(passkey, salt=b"PyRosettaCluster")

    return hkdf_expand(prk, info)
