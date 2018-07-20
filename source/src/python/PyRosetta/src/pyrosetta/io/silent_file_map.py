import os
from os import path
import warnings


def _open_gzip_compressed(f):
    from gzip import GzipFile
    return GzipFile(filename=f) if isinstance(f, str) else GzipFile(fileobj=f)

def _open_bz2_compressed(f):
    from bz2file import BZ2File
    return BZ2File(f)

# bz2file used to support concatenated bz2 streams
_compressed_magic_bytes = { }
_compressed_magic_bytes[b"\x1f\x8b\x08"] = _open_gzip_compressed
_compressed_magic_bytes[b"\x42\x5a\x68"] = _open_bz2_compressed


def open_compressed(filename):
    """Open a possibly compressed file or stream as a decompressed file object.
    Checks for gzip/bz2 prefix bytes and returns a decompression stream if needed.
    filename - str filename or file object supporting seek
    """

    if isinstance(filename, str):
        file_object = open(filename, "rb")
    else:
        file_object = filename

    # Peek ahead for magic bytes
    current_stream_position = file_object.tell()
    file_prefix = file_object.read(max(len(b) for b in _compressed_magic_bytes))
    file_object.seek(current_stream_position)

    for magic_bytes in _compressed_magic_bytes:
        if file_prefix.startswith(magic_bytes):
            return _compressed_magic_bytes[magic_bytes](file_object)
    else:
        return file_object

from io import BytesIO
try:
    from collections.abc import MutableMapping, Mapping
except ImportError:
    # these types are in the collections module in python < 3.3
    from collections import MutableMapping, Mapping

class SilentFilePoseAccessor(Mapping):
    """Map class providing read-only pose access to protein silent file data."""

    def __init__(self, sdata):
        """Initialize pose accessor
        Args:
            sdata - SilentFileData instance
        """
        self.sdata = sdata

    def keys(self):
        """Return list of all tags in silent file data."""
        return list(self.sdata.tags())

    def __getitem__(self, tag_name):
        """Get pose corresponding to given tag."""
        if not self.sdata.has_tag(tag_name):
            raise KeyError(tag_name)

        import pyrosetta.rosetta.core.pose

        p = pyrosetta.rosetta.core.pose.Pose()
        sfd = self.sdata.get_structure(tag_name)

        sfd.fill_pose(p)
        p.pdb_info().name( tag_name )

        return p

    def __len__(self):
        return self.keys().__len__()

    def __iter__(self):
        return self.keys().__iter__()

class SilentFileEnergyAccessor(Mapping):
    """Map class providing read-only energies access to protein silent file data."""

    def __init__(self, sdata):
        """Initialize energies accessor
        Args:
            sdata - SilentFileData instance
        """
        self.sdata = sdata

    def keys(self):
        """Return list of all tags in silent file data."""
        return list(self.sdata.tags())

    def __getitem__(self, tag_name):
        """Get energies corresponding to given tag."""
        if not self.sdata.has_tag(tag_name):
            raise KeyError(tag_name)

        sfd = self.sdata.get_structure(tag_name)
        return { e : sfd.get_energy(e) for e in sfd.energy_names().energy_names() }

    def __len__(self):
        return self.keys().__len__()

    def __iter__(self):
        return self.keys().__iter__()

class SilentFileMap(MutableMapping):
    """Map class providing key-value interface to silent file data."""

    def __init__(self, filename = None, data = None):
        """Initialize mapping over given filename or silent file data.
        filename - Read silent data from given file if data not provided, default path for `write`.
        data - String silent file data to load.
        """
        import pyrosetta
        import pyrosetta.rosetta.core.io.silent as silent

        self.filename = filename
        self.file_writable = True

        if data is None and self.filename is not None and path.lexists(self.filename):
            with open_compressed(filename) as opened_file:
                if not isinstance(opened_file, BytesIO):
                    # Can not rewrite compressed files
                    # open_compressed returns 'file' object for uncompressed input streams
                    file_writable = False
                data = opened_file.read()

        if data:
            sdata = silent.SilentFileData((self.filename if self.filename else "<anonymous>"), silent.SilentFileOptions())
            # Need iostring wrapper?
            sdata.read_stream(
                pyrosetta.rosetta.std.istringstream(data),
                pyrosetta.rosetta.utility.vector1_string(),
                False)
        else:
            sdata = silent.SilentFileData((self.filename if self.filename else "<anonymous>"), silent.SilentFileOptions())

        self.sdata = sdata

    @property
    def tags(self):
        """Return list of all tags in silent file data."""
        return list(self.sdata.tags())

    @property
    def poses(self):
        """Return SilentFilePoseAccessor over silent file map."""
        return SilentFilePoseAccessor(self.sdata)

    @property
    def energies(self):
        """Return SilentFileEnergiesAccessor over silent file map."""
        return SilentFileEnergyAccessor(self.sdata)

    def write_file(self, filename=None, force=False):
        """Write file to given filename."""

        if filename is None:
            if self.filename and not self.file_writable:
                raise ValueError("Unable to write target file: %r" % self.filename)
            elif not self.filename:
                raise ValueError("Unable to resolve target file.")
            else:
                filename = self.filename

        if path.lexists(filename):
            if not force:
                raise ValueError("Target filename exists: %s" % filename)
            else:
                os.remove(filename)

        self.sdata.write_all(filename)

    def __getitem__(self, tag_name):
        if not self.sdata.has_tag(tag_name):
            raise KeyError(tag_name)

        return self.sdata.get_structure(tag_name)

    def __setitem__(self, tag_name, structure):
        import pyrosetta.rosetta.core.pose as pose
        import pyrosetta.rosetta.core.io.silent as silent

        if not isinstance(tag_name, str):
            raise ValueError(tag_name)

        if self.sdata.has_tag(tag_name):
            raise NotImplementedError("__setitem__ not implemented for existing keys")

        if isinstance(structure, pose.Pose):
            structure = silent.BinarySilentStruct(silent.SilentFileOptions(), structure, tag_name)
        elif isinstance(structure, silent.SilentStruct):
            structure = structure.clone()
            structure.decoy_tag(tag_name)
        else:
            raise ValueError("Unable to store structure type: %s" % structure)

        self.sdata.add_structure(structure)

    def __delitem__(self, tag_name):
        raise NotImplementedError("__delitem__ not implemented.")

    def __len__(self):
        return self.tags.__len__()

    def __iter__(self):
        return self.tags.__iter__()

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        if not type:
            self.write_file(force=True)
