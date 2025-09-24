# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import json
import hmac
import pyrosetta
import pyrosetta.distributed.io as io

from pyrosetta.distributed.packed_pose.core import PackedPose
from typing import (
    Any,
    Dict,
    Generic,
    List,
    Optional,
    Tuple,
    TypeVar,
)

from pyrosetta.distributed.cluster.hkdf import HASHMOD, compare_digest, derive_init_key


G = TypeVar("G")


class InitFileSigner(Generic[G]):
    """Sign or verify PyRosetta initialization files by PyRosettaCluster."""

    _encoding = "utf-8"
    _prefix = b'PyRosettaCluster_init_file_signer'

    def __init__(self, input_packed_pose=None, output_packed_pose=None, metadata=None) -> None:
        self.inp_pkl = self._to_pickle(input_packed_pose)
        self.out_pkl = self._to_pickle(output_packed_pose)
        self.metadata_pkl = self._to_encoding(metadata)

    def _to_pickle(self, packed_pose: Optional[PackedPose]) -> bytes:
        return io.to_pickle(packed_pose) if isinstance(packed_pose, PackedPose) else (b'\x00' * 32)

    def _to_encoding(self, obj: Any) -> bytes:
        return json.dumps(
            obj,
            skipkeys=False,
            ensure_ascii=False,
            check_circular=True,
            allow_nan=False,
            cls=None,
            indent=None,
            separators=(",", ":"),
            default=None,
            sort_keys=True, # Deterministic
        ).encode(InitFileSigner._encoding)

    def _get_pose_digest(self, pkl: bytes) -> bytes:
        return HASHMOD(pkl).digest()

    def _join_bytes(self, *values: List[bytes]) -> bytes:
        return b'+'.join(values)

    def _setup_poses_pair(self, inp_pkl: bytes, out_pkl: bytes) -> bytes:
        return self._join_bytes(InitFileSigner._prefix, inp_pkl, out_pkl)

    def _get_poses_digest(self, inp_pkl: bytes, out_pkl: bytes) -> bytes:
        return HASHMOD(self._setup_poses_pair(inp_pkl, out_pkl)).digest()

    def _get_poses_hexdigest(self, inp_pkl: bytes, out_pkl: bytes) -> str:
        return HASHMOD(self._setup_poses_pair(inp_pkl, out_pkl)).hexdigest()

    def _get_pkg_data(self) -> bytes:
        return self._join_bytes(
            InitFileSigner._prefix,
            pyrosetta._version_string().encode(InitFileSigner._encoding),
            pyrosetta.distributed.cluster.__version__.encode(InitFileSigner._encoding),
        )

    def _get_hmac_hexdigest(self, key: bytes, data: bytes) -> str:
        return hmac.new(key, data, HASHMOD).hexdigest()

    def _get_init_key_and_msg(self) -> Tuple[bytes, bytes]:
        inp_digest = self._get_pose_digest(self.inp_pkl)
        out_digest = self._get_pose_digest(self.out_pkl)
        poses_digest = self._get_poses_digest(self.inp_pkl, self.out_pkl)
        pkg_data = self._join_bytes(self._get_pkg_data(), self.metadata_pkl, poses_digest)
        init_key = derive_init_key(b'InitFileSigner', pkg_data)
        msg = self._join_bytes(InitFileSigner._prefix, inp_digest, out_digest, pkg_data)

        return init_key, msg

    def sign_sha256(self) -> str:
        """Sign PyRosetta initialization file pose data."""
        return self._get_poses_hexdigest(self.inp_pkl, self.out_pkl)

    def sign_digest(self) -> str:
        """Sign PyRosetta initialization file package data, pose data, and metadata."""
        return self._get_hmac_hexdigest(*self._get_init_key_and_msg())

    def sign(self) -> Dict[str, str]:
        """Return a `dict` object with SHA256 and HMAC signature metadata."""
        return {
            "sha256": self.sign_sha256(),
            "signature": self.sign_digest(),
        }

    def verify_sha256(self, sha256: Optional[str]) -> bool:
        """Verify PyRosetta initialization file pose data."""
        return sha256 == self._get_poses_hexdigest(self.inp_pkl, self.out_pkl)

    def verify_signature(self, signature: Optional[str]) -> bool:
        """Verify PyRosetta initialization file package data, pose data, and metadata."""
        return compare_digest(signature, self.sign_digest())

    def verify(self, sha256: Optional[str], signature: Optional[str]) -> bool:
        """Verify PyRosetta initialization file SHA256 and HMAC signature metadata."""
        return self.verify_sha256(sha256) and self.verify_signature(signature)
