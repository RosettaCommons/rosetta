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
import struct

from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.rosetta.core.pose import (
    Pose,
    get_all_comments,
)
from typing import (
    Any,
    Dict,
    Generic,
    List,
    NoReturn,
    Optional,
    Tuple,
    TypeVar,
    Union,
)

from pyrosetta.distributed.cluster.hkdf import (
    HASHMOD,
    compare_digest,
    derive_init_key,
)

G = TypeVar("G")


class PackedPoseHasher(Generic[G]):
    """Digest the scientific state of a `PackedPose` or `Pose` object."""

    _encoding: str = "utf-8"
    _default_bytes: bytes = (b'\x00' * 32)

    def __init__(
        self,
        packed_pose: Optional[Union[PackedPose, Pose]] = None,
        include_cache: bool = False,
        include_comments: bool = False,
    ) -> None:
        """
        Initialize the `PackedPoseHasher` class.

        Args:
            `packed_pose`:
                A `PackedPose` or `Pose` object to hash. If `None`, then other keyword arguments have no effect
                and a constant hash is returned from the `PackedPoseHasher.digest` method.

                Default: `None`

            `include_cache`:
                A `bool` object specifying whether or not to include the `Pose.cache` entries.

                Default: `False`

             `include_comments`:
                A `bool` object specifying whether or not to include the `Pose` comments entries.

                Default: `False`

        Returns:
            `None`
        """
        self.pose = io.to_pose(packed_pose)
        self.include_cache = include_cache
        self.include_comments = include_comments
        self.hashmod = HASHMOD()

    def digest(self) -> bytes:
        """Digest the `PackedPose` or `Pose` object, otherwise return a default value."""
        if isinstance(self.pose, Pose):
            self.add_coordinates()
            if self.include_cache:
                self.add_cache()
            if self.include_comments:
                self.add_comments()
            return self.hashmod.digest()
        else:
            return PackedPoseHasher._default_bytes

    def encode_string(self, obj: str) -> bytes:
        """Encode a string object into bytes."""
        return obj.strip().encode(PackedPoseHasher._encoding)

    def update_hashmod(self, obj: Any) -> Optional[NoReturn]:
        """Update the hashmod with an input object."""
        if isinstance(obj, str):
            self.hashmod.update(self.encode_string(obj))
        elif isinstance(obj, int):
            self.hashmod.update(struct.pack("!i", obj))
        elif isinstance(obj, float):
            self.hashmod.update(struct.pack("!d", obj))
        elif isinstance(obj, bytes):
            self.hashmod.update(obj)
        elif isinstance(obj, bytearray):
            self.hashmod.update(bytes(obj))
        else:
            raise NotImplementedError(f"Object type not supported: {type(obj)}")

    def add_coordinates(self) -> None:
        """
        Update hashmod with residue numbers, residue names, atom numbers, atom names, and double precision
        atomic coordinate components of the `Pose` object.
        """
        for res in range(1, self.pose.size() + 1):
            residue = self.pose.residue(res)
            # Add residue number
            self.update_hashmod(res)
            # Add residue name
            self.update_hashmod(residue.name())
            for atom in range(1, residue.natoms() + 1):
                # Add atom number
                self.update_hashmod(atom)
                # Add atom name
                self.update_hashmod(residue.atom_name(atom))
                # Add atom coordinate components
                xyz = residue.atom(atom).xyz()
                for axis in "xyz":
                    self.update_hashmod(getattr(xyz, axis))

    def add_cache(self) -> None:
        """Update hashmod with serialized `Pose.cache` dictionary entries."""
        for entry in (
            self.pose.cache.energies.all,
            self.pose.cache.extra.string.all,
            self.pose.cache.extra.real.all,
            self.pose.cache.metrics.string.all,
            self.pose.cache.metrics.real.all,
            self.pose.cache.metrics.composite_string.all,
            self.pose.cache.metrics.composite_real.all,
            self.pose.cache.metrics.per_residue_string.all,
            self.pose.cache.metrics.per_residue_real.all,
        ):
            for key, value in sorted(entry.items()):
                self.update_hashmod(key)
                self.update_hashmod(value)
        for key, value in sorted(self.pose.cache.metrics.per_residue_probabilities.all.items()):
            self.update_hashmod(key)
            for name3, prob in sorted(value.items()):
                self.update_hashmod(name3)
                self.update_hashmod(prob)

    def add_comments(self) -> None:
        """Update hashmod with raw `Pose` comments."""
        comments = dict(get_all_comments(self.pose))
        for key, value in sorted(comments.items()):
            self.update_hashmod(key)
            self.update_hashmod(value)


class InitFileSigner(Generic[G]):
    """Sign or verify PyRosetta initialization files by `PyRosettaCluster`."""

    _encoding: str = "utf-8"
    _prefix: bytes = b'PyRosettaCluster_init_file_signer'

    def __init__(self, input_packed_pose=None, output_packed_pose=None, metadata=None) -> None:
        self.inp_pkl = self._to_packed_pose_hash(input_packed_pose)
        self.out_pkl = self._to_packed_pose_hash(output_packed_pose)
        self.metadata_pkl = self._to_encoding(metadata)

    def _to_packed_pose_hash(self, packed_pose: Optional[PackedPose]) -> bytes:
        return PackedPoseHasher(packed_pose, include_cache=True, include_comments=True).digest()

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
