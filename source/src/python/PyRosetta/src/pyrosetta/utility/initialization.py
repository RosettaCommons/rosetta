from __future__ import absolute_import
# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

###############################################################################
# PyRosetta initialization files


__author__ = "Jason C. Klima"


import base64
import collections
import datetime
import hashlib
import inspect
import json
import os
import pyrosetta
import re
import tempfile
import warnings
import zlib

from pprint import pprint
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.rosetta.core.simple_metrics.composite_metrics import ProtocolSettingsMetric


class PyRosettaInitFileParserBase(object):
    _database_option_name = "in:path:database"
    _corrections_option_name = "corrections:"
    _init_file_extension = ".init"
    _strftime_format = "%Y-%m-%d-%H-%M-%S"

    def get_pyrosetta_build(self):
        return pyrosetta._version_string()

    def pyrosetta_build_warning(self, original_pyrosetta_build, current_pyrosetta_build):
        _msg = os.linesep.join(
            [
                "The PyRosetta version that generated the '{0}' file ".format(self._init_file_extension)
                + "does not match the current PyRosetta version. Please inspect "
                + "the PyRosetta input files if you encounter any issues during "
                + "or after PyRosetta initialization: {0}".format(self.kwargs["output_dir"]),
                "Original: {0}".format(original_pyrosetta_build),
                "Current:  {0}".format(current_pyrosetta_build),
            ]
        )
        warnings.warn(_msg, UserWarning, stacklevel=4)

    @property
    def has_cereal(self):
        try:
            import pyrosetta.rosetta.cereal as cereal  # noqa: F401
            assert hasattr(cereal, "BinaryInputArchive")
            assert hasattr(cereal, "BinaryOutputArchive")
            return True
        except (ImportError, AssertionError):
            return False

    def get_to_base64(self):
        assert self.has_cereal, (
            f"To cache `Pose` and/or `PackedPose` objects in the output '{self._init_file_extension}' "
            + "file, please ensure that PyRosetta is built with serialization support. "
            + "Current PyRosetta build: {0}".format(self.get_pyrosetta_build())
        )
        try:
            import pyrosetta.distributed.io as io  # noqa: F401
            return io.to_base64
        except ImportError as ex:
            raise ImportError(
                "{0}. Please ensure that all 'pyrosetta.distributed' framework ".format(ex)
                + "package dependencies are installed. For instructions, please visit:\n"
                + "https://pypi.org/project/pyrosetta-installer\n"
                + "https://pypi.org/project/pyrosetta-distributed\n"
            )

    @property
    def was_init_called(self):
        return pyrosetta.rosetta.basic.was_init_called()

    def validate_init_was_called(self):
        if not self.was_init_called:
            raise RuntimeError(
                "PyRosetta must be already initialized to dump a '{0}' file. ".format(self._init_file_extension)
                + "Please run `pyrosetta.init()` with custom options and try again."
            )

    def validate_init_was_not_called(self):
        if self.was_init_called:
            raise RuntimeError(
                "PyRosetta must not be already initialized to initialize from a '{0}' file. ".format(self._init_file_extension)
                + "Please ensure that `pyrosetta.init()` was not already called (e.g., "
                + "if using a Jupyter notebook, please restart the kernel) and try again."
            )

    def md5_warning(self, md5, expected_md5):
        if md5 is None:
            warnings.warn(
                "The 'md5' key is missing from the '{0}' file!".format(self._init_file_extension),
                UserWarning,
                stacklevel=5,
            )
        elif isinstance(md5, str):
            if md5 != expected_md5:
                warnings.warn(
                    (
                        "The expected MD5 value differs from the 'md5' key value in the '{0}' file!\n".format(self._init_file_extension)
                        + "Expected: '{0}'\n".format(expected_md5)
                        + "Value:    '{0}'\n".format(md5)
                    ),
                    UserWarning,
                    stacklevel=5,
                )
        else:
            raise TypeError(self._malformed_init_file_error_msg)


class PyRosettaInitFileSerializer(object):
    _chunk_size = 64 * 1024 # bytes
    _compression_level = 9
    _encoding = "utf-8"
    _prefix_string = "[PyRosettaInitTextFile]"
    _prefix_binary = "[PyRosettaInitBinaryFile]"
    _tag_str = b'Txt'
    _tag_obj = b'Obj'

    def dump_json(self, obj):
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
            sort_keys=False,
        )

    def load_json(self, string):
        return json.loads(
            string,
            cls=None,
            object_hook=None,
            object_pairs_hook=None,
            parse_int=None,
            parse_constant=None,
        )

    def join_tag(self, tag, raw):
        return tag + raw

    def split_tag(self, obj):
        for t in (PyRosettaInitFileSerializer._tag_str, PyRosettaInitFileSerializer._tag_obj):
            if obj.startswith(t):
                tag, raw = obj[:len(t)], obj[len(t):]
                break
        else:
            ValueError(obj)

        return tag, raw

    def encode_bytestring(self, bytestring):
        return base64.b64encode(bytestring).decode(PyRosettaInitFileSerializer._encoding, errors="strict")

    def encode_string(self, obj):
        if isinstance(obj, str):
            tag = PyRosettaInitFileSerializer._tag_str
            raw = obj
        elif isinstance(obj, (list, dict)):
            tag = PyRosettaInitFileSerializer._tag_obj
            raw = self.dump_json(obj)
        else:
            raise TypeError(obj)
        compressed = zlib.compress(
            raw.encode(PyRosettaInitFileSerializer._encoding),
            PyRosettaInitFileSerializer._compression_level,
        )

        return self.encode_bytestring(self.join_tag(tag, compressed))

    def decode_binary(self, string):
        return base64.b64decode(string, validate=True)

    def decode_string(self, bytestring, max_decompressed_bytes):
        obj = self.decode_binary(bytestring)
        tag, raw = self.split_tag(obj)
        decompressed = self.zlib_decompress(raw, max_decompressed_bytes).decode(
            PyRosettaInitFileSerializer._encoding, errors="strict"
        )
        if tag == PyRosettaInitFileSerializer._tag_str:
            result = decompressed
        elif tag == PyRosettaInitFileSerializer._tag_obj:
            result = self.load_json(decompressed)
        else:
            raise ValueError(tag)

        return result

    def zlib_decompress(self, data, max_decompressed_bytes):
        buf = memoryview(data)
        zobj = zlib.decompressobj()
        arr = bytearray()
        for i in range(0, len(buf), PyRosettaInitFileSerializer._chunk_size):
            arr += zobj.decompress(buf[i: (i + PyRosettaInitFileSerializer._chunk_size)])
            if len(arr) > max_decompressed_bytes:
                raise BufferError(self.get_zlib_decompress_err_msg(arr, max_decompressed_bytes))
        arr += zobj.flush()
        if len(arr) > max_decompressed_bytes:
            raise BufferError(self.get_zlib_decompress_err_msg(arr, max_decompressed_bytes))

        return bytes(arr)

    def get_zlib_decompress_err_msg(self, arr, max_decompressed_bytes):
        return "Decompressed data exceeds maximum bytes size limit: {0} > {1}".format(len(arr), max_decompressed_bytes)

    @staticmethod
    def get_md5(init_dict):
        return hashlib.md5(
            json.dumps(
                init_dict,
                skipkeys=False,
                ensure_ascii=False,
                check_circular=True,
                allow_nan=False,
                cls=None,
                indent=None,
                separators=(",", ":"),
                default=None,
                sort_keys=True,
            ).encode(PyRosettaInitFileSerializer._encoding)
        ).hexdigest()


class PyRosettaInitFileWriter(PyRosettaInitFileParserBase, PyRosettaInitFileSerializer):
    def __init__(self, output_filename, **kwargs):
        self.validate_init_was_called()
        self.kwargs = self.setup_kwargs(**kwargs)
        self.output_filename = self.setup_output_filename(output_filename)
        self.cached_files = []

    def setup_output_filename(self, output_filename):
        if not isinstance(output_filename, str):
            raise TypeError(
                "Output file must be a `str` object. Received: {0}".format(type(output_filename))
            )
        if not output_filename.endswith(self._init_file_extension):
            raise NameError(
                "Output file must end with the '{0}' filename extension.".format(self._init_file_extension)
            )
        if os.path.isfile(output_filename) and not self.kwargs["overwrite"]:
            raise FileExistsError(
                "Output '{0}' file already exists! Please remove the file and try again: {1}".format(
                    self._init_file_extension, output_filename
                )
            )

        return output_filename

    def setup_kwargs(self, **kwargs):
        if "poses" in kwargs and kwargs["poses"] is None:
            kwargs["poses"] = []
        else:
            to_base64 = self.get_to_base64()
            objs = kwargs["poses"]
            if isinstance(objs, (Pose, PackedPose)):
                kwargs["poses"] = [to_base64(objs)]
            elif isinstance(objs, collections.abc.Iterable):
                kwargs["poses"] = []
                for obj in objs:
                    if isinstance(obj, (Pose, PackedPose)):
                        kwargs["poses"].append(to_base64(obj))
                    else:
                        raise TypeError(
                            "The 'poses' keyword argument parameter must be a `Pose` or `PackedPose` object, "
                            + "or an iterable of `Pose` or `PackedPose` objects. Received: {0}".format(type(obj))
                        )
        for key in ("author", "email", "license"):
            if key in kwargs and kwargs[key] is None:
                kwargs[key] = ""
            elif not isinstance(kwargs[key], str):
                raise TypeError(
                    "The '{0}' keyword argument parameter must be a `str` object. Received: {1}".format(
                        key, type(kwargs[key])
                    )
                )
        if "metadata" in kwargs and kwargs["metadata"] is None:
            kwargs["metadata"] = {}
        else:
            self.assert_metadata_json_serializable(kwargs["metadata"])
        kwargs["pyrosetta_build"] = self.get_pyrosetta_build()
        kwargs["datetime"] = self.get_datetime_now()
        if "overwrite" in kwargs and kwargs["overwrite"] is None:
            kwargs["overwrite"] = False
        elif not isinstance(kwargs["overwrite"], bool):
            raise TypeError(
                "The 'overwrite' keyword argument parameter must be a `bool` object. Received: {0}".format(type(kwargs["overwrite"]))
            )
        if "dry_run" in kwargs and kwargs["dry_run"] is None:
            kwargs["dry_run"] = False
        elif not isinstance(kwargs["dry_run"], bool):
            raise TypeError(
                "The 'dry_run' keyword argument parameter must be a `bool` object. Received: {0}".format(type(kwargs["dry_run"]))
            )
        if "verbose" in kwargs and kwargs["verbose"] is None:
            kwargs["verbose"] = True
        elif not isinstance(kwargs["verbose"], bool):
            raise TypeError(
                "The 'verbose' keyword argument parameter must be a `bool` object. Received: {0}".format(type(kwargs["verbose"]))
            )

        return kwargs

    def assert_metadata_json_serializable(self, data):
        try:
            json.dumps(data)
        except:
            raise ValueError("Input 'metadata' keyword argument parameter must be JSON-serializable.")

    def get_datetime_now(self):
        return datetime.datetime.now(datetime.timezone.utc).strftime(self._strftime_format)

    def init_pose(self):
        return Pose()

    def get_protocol_settings_metric(
        self,
        base_name_option_only=False,
        get_script_vars=False,
        get_user_options=True,
        skip_corrections=False,
    ):
        metric = ProtocolSettingsMetric()
        options = pyrosetta.rosetta.basic.options.process()
        metric.parse_options(
            options=options,
            base_name_option_only=base_name_option_only,
            get_script_vars=get_script_vars,
            get_user_options=get_user_options,
            skip_corrections=skip_corrections,
        )

        return metric

    def get_protocol_settings_dict(self):
        pose = self.init_pose()
        metric = self.get_protocol_settings_metric()

        return dict(metric.calculate(pose))

    def get_options_dict(self):
        options_dict = {}
        for option_name, string in self.get_protocol_settings_dict().items():
            values = string.split()
            if option_name == self._database_option_name:
                options_dict[option_name] = sorted(set(values), key=values.index)
            else:
                options_dict[option_name] = values

        return options_dict

    def get_options_str(self):
        options_dict = self.get_options_dict()

        return " ".join(
            [
                "-{0} {1}".format(option_name, " ".join(values))
                for option_name, values in options_dict.items()
            ]
        )

    def get_encoded_options_dict(self):
        options_dict = self.get_options_dict()
        encoded_options_dict = collections.defaultdict(list)
        for option_name, values in options_dict.items():
            for value in values:
                if os.path.isfile(value):
                    encoded_options_dict[option_name].append(self.encode_file(value))
                elif os.path.isdir(value):
                    rel_value = os.path.relpath(value, start=os.curdir) if os.path.isabs(value) else value
                    if value != rel_value and option_name != self._database_option_name:
                        warnings.warn(
                            "The option '-{0}' with path '{1}' is being cached as the relative path: '{2}'.".format(
                                option_name, value, rel_value
                            ),
                            UserWarning,
                            stacklevel=1,
                        )
                    encoded_options_dict[option_name].append(rel_value)
                else:
                    encoded_options_dict[option_name].append(value)

        return dict(encoded_options_dict)

    def is_text_file(self, filename):
        try:
            with open(filename, "r") as f:
                for _line in f:
                    return True
        except UnicodeDecodeError:
            return False

    def encode_object(self, obj):
        return self.encode_string(obj)

    def format_encode_bytestring(self, bytestring):
        return "{0}{1}".format(
            PyRosettaInitFileSerializer._prefix_binary,
            self.encode_bytestring(bytestring),
        )

    def format_encode_string(self, string, parent_dir):
        obj = self.encode(string, parent_dir)

        return "{0}{1}".format(
            PyRosettaInitFileSerializer._prefix_string,
            self.encode_object(obj),
        )

    def format_encode_substring(self, string):
        return "{0}{1}".format(
            PyRosettaInitFileSerializer._prefix_string,
            self.encode_string(string),
        )

    def encode_subfile(self, file):
        basename = os.path.basename(file)
        if self.is_text_file(file):
            with open(file, "r") as f:
                result = {basename: self.format_encode_substring(f.read())}
        else:
            with open(file, "rb") as f:
                result = {basename: self.format_encode_bytestring(f.read())}

        return result

    def encode(self, string, parent_dir):
        results = []
        for obj in filter(bool, re.split(r'(\s+)', string)):
            supported_subfiles = {
                obj,  # Full path to a subfile
                os.path.join(parent_dir, obj),  # Subfile may be in the directory of the parent file
            }
            for file in supported_subfiles:
                if os.path.isfile(file):
                    self.cached_files.append(file)
                    result: dict = self.encode_subfile(file)  # Reserve `dict` object for a subfile
                    break
            else:  # Not a subfile
                result: str = obj  # Reserve `str` object for text
            results.append(result)

        return results

    def encode_file(self, filename):
        self.cached_files.append(filename)
        results = {}
        if self.is_text_file(filename):
            parent_dir = os.path.dirname(filename)
            with open(filename, "r") as f:
                result = self.format_encode_string(f.read(), parent_dir)
        else:
            with open(filename, "rb") as f:
                result = self.format_encode_bytestring(f.read())
        results[os.path.basename(filename)] = result

        return results

    def print_cached_files(self, dry_run):
        if dry_run:
            print(f"Dry run dump PyRosetta '{self._init_file_extension}' file:")
        else:
            print(f"Dumping PyRosetta '{self._init_file_extension}' file to: {self.output_filename}")
        if len(self.kwargs["poses"]) > 0:
            print("Compressed {0} PyRosetta `PackedPose` object(s).".format(len(self.kwargs["poses"])))
        else:
            print("No input PyRosetta `Pose` or `PackedPose` object(s) to compress.")
        if len(self.cached_files) > 0:
            print("Compressed {0} PyRosetta initialization input file(s):".format(len(self.cached_files)))
            for file in self.cached_files:
                print(os.path.relpath(file, start=os.curdir))
        else:
            print("No PyRosetta initialization input files to compress.")
        if dry_run:
            print(f"Skipping dumping PyRosetta initialization '{self._init_file_extension}' file...")

    @staticmethod
    def write_json(data_dict, output_filename):
        with open(output_filename, "w") as f:
            return json.dump(
                data_dict,
                f,
                skipkeys=False,
                ensure_ascii=False,
                check_circular=True,
                allow_nan=False,
                cls=None,
                indent=2,
                separators=(", ", ": "),
                default=None,
                sort_keys=True,
            )

    def dump(self):
        encoded_options_dict = self.get_encoded_options_dict()
        encoded_options_dict.pop(self._database_option_name, None)
        overwrite = self.kwargs.pop("overwrite")
        dry_run = self.kwargs.pop("dry_run")
        verbose = self.kwargs.pop("verbose")
        self.kwargs["num_files"] = len(self.cached_files)
        self.kwargs["num_poses"] = len(self.kwargs["poses"])
        data_dict = {
            **self.kwargs,
            "options": encoded_options_dict,
        }
        data_dict["md5"] = PyRosettaInitFileSerializer.get_md5(data_dict)
        if verbose:
            self.print_cached_files(dry_run)
        if not dry_run and ((not os.path.isfile(self.output_filename)) or overwrite):
            PyRosettaInitFileWriter.write_json(data_dict, self.output_filename)
            if verbose:
                print(
                    f"Dumped PyRosetta '{self._init_file_extension}' file size:",
                    round(os.path.getsize(self.output_filename) * 1e-6, 3),
                    "MB",
                    sep=" ",
                )


class PyRosettaInitFileReader(PyRosettaInitFileParserBase, PyRosettaInitFileSerializer):
    def __init__(self, init_file, **kwargs):
        self.init_file = init_file
        self.init_dict = self.setup_init_dict(init_file)
        self.kwargs = self.setup_kwargs(**kwargs)
        self.file_counter = 0

    def setup_kwargs(self, **kwargs):
        if kwargs["dry_run"] is None:
            kwargs["dry_run"] = False
        if not isinstance(kwargs["dry_run"], bool):
            raise TypeError("The 'dry_run' keyword argument parameter must be a `bool` object. Received: {0}".format(type(kwargs["dry_run"])))
        if kwargs["output_dir"] is None:
            output_dir = os.path.join(os.getcwd(), "pyrosetta_init_files")
        elif isinstance(kwargs["output_dir"], str):
            output_dir = os.path.abspath(kwargs["output_dir"])
        else:
            raise TypeError("The 'output_dir' keyword argument parameter must be a `str` object. Received: {0}".format(type(kwargs["output_dir"])))
        if not os.path.isdir(output_dir):
            kwargs["output_dir"] = output_dir
        elif kwargs["dry_run"]:
            warnings.warn(
                "The output directory already exists! Please remove the output directory before disabling dry run: {0}".format(output_dir),
                UserWarning,
                stacklevel=4,
            )
        else:
            raise IsADirectoryError(
                "The output directory already exists! Please remove the output directory and try again: {0}".format(output_dir)
            )
        if not isinstance(kwargs["skip_corrections"], (bool, type(None))):
            raise TypeError(
                "The 'skip_corrections' keyword argument parameter must be a `bool` or `NoneType` object. Received: {0}".format(type(kwargs["skip_corrections"]))
            )
        if kwargs["relative_paths"] is None:
            kwargs["relative_paths"] = False
        if not isinstance(kwargs["relative_paths"], bool):
            raise TypeError(
                "The 'relative_paths' keyword argument parameter must be a `bool` object. Received: {0}".format(type(kwargs["relative_paths"]))
            )
        if kwargs["max_decompressed_bytes"] is None:
            kwargs["max_decompressed_bytes"] = 200_000_000 # 200 MB
        if not isinstance(kwargs["max_decompressed_bytes"], int):
            raise TypeError(
                "The 'max_decompressed_bytes' keyword argument parameter must be a `int` object. Received: {0}".format(type(kwargs["max_decompressed_bytes"]))
            )
        elif kwargs["max_decompressed_bytes"] <= 0:
            raise ValueError(
                "The 'max_decompressed_bytes' keyword argument parameter must be greater than 0 bytes. Received: {0}".format(kwargs["max_decompressed_bytes"])
            )
        if kwargs["database"] is None:
            kwargs["database"] = pyrosetta._rosetta_database_from_env()
        if not (isinstance(kwargs["database"], str) and os.path.isdir(kwargs["database"])):
            raise NotADirectoryError("PyRosetta database directory not found: {0}".format(kwargs["database"]))
        if kwargs["verbose"] is None:
            kwargs["verbose"] = True
        if not isinstance(kwargs["verbose"], bool):
            raise TypeError(
                "The 'verbose' keyword argument parameter must be a `bool` object. Received: {0}".format(type(kwargs["verbose"]))
            )
        fullargspec = inspect.getfullargspec(pyrosetta.init)
        default_init_kwargs = dict(zip(fullargspec.args, fullargspec.defaults))
        if kwargs["set_logging_handler"] is None:
            kwargs["set_logging_handler"] = default_init_kwargs["set_logging_handler"]
        if kwargs["notebook"] is None:
            kwargs["notebook"] = default_init_kwargs["notebook"]
        if kwargs["silent"] is None:
            kwargs["silent"] = default_init_kwargs["silent"]

        return kwargs

    @staticmethod
    def read_json(init_file):
        with open(init_file, "r") as f:
            return json.load(
                f,
                cls=None,
                object_hook=None,
                parse_float=None,
                parse_int=None,
                parse_constant=None,
                object_pairs_hook=None,
            )

    def setup_init_dict(self, init_file):
        if isinstance(init_file, str) and init_file.endswith(self._init_file_extension) and os.path.isfile(init_file):
            _init_dict = PyRosettaInitFileReader.read_json(init_file)
            _md5 = _init_dict.pop("md5", None)
            _expected_md5 = PyRosettaInitFileSerializer.get_md5(_init_dict)
            self.md5_warning(_md5, _expected_md5)
            _init_dict.pop("poses", None)
            return _init_dict
        else:
            raise ValueError(
                "Please provide a valid input PyRosetta initialization '{0}' file. Received: {1}".format(self._init_file_extension, init_file)
            )

    def get_encoded_options_dict(self):
        encoded_options_dict = self.init_dict["options"]
        encoded_options_dict[self._database_option_name] = [self.kwargs["database"]]

        return encoded_options_dict

    @property
    def _malformed_init_file_error_msg(self):
        return "Cannot read malformed PyRosetta '{0}' file: {1}".format(self._init_file_extension, self.init_file)

    def format_decode_binary(self, value):
        return self.decode_binary(value.split(PyRosettaInitFileSerializer._prefix_binary)[-1])

    def format_decode_string(self, value, option_name):
        obj = self.decode_string(value.split(PyRosettaInitFileSerializer._prefix_string)[-1], self.kwargs["max_decompressed_bytes"])

        return self.decode(obj, option_name)

    def format_decode_substring(self, value):
        return self.decode_string(value.split(PyRosettaInitFileSerializer._prefix_string)[-1], self.kwargs["max_decompressed_bytes"])

    def decode(self, encoded_object, option_name):
        file_content = ""
        for obj in encoded_object:
            if isinstance(obj, dict):  # a `dict` object is reserved for a subfile
                assert len(obj) == 1, self._malformed_init_file_error_msg
                basename, data = next(iter(obj.items()))
                assert data.startswith(
                    (PyRosettaInitFileSerializer._prefix_string, PyRosettaInitFileSerializer._prefix_binary)
                ), self._malformed_init_file_error_msg
                if data.startswith(PyRosettaInitFileSerializer._prefix_string):
                    subfile_content = self.format_decode_substring(data)
                    subfilename = self.write_text_file(option_name, basename, subfile_content)
                elif data.startswith(PyRosettaInitFileSerializer._prefix_binary):
                    subfile_content = self.format_decode_binary(data)
                    subfilename = self.write_binary_file(option_name, basename, subfile_content)
                result: str = subfilename
            elif isinstance(obj, str):  # a `str` object is reserved for text
                assert bool(obj), self._malformed_init_file_error_msg
                result: str = obj
            else:
                raise ValueError(self._malformed_init_file_error_msg)
            file_content += result

        return file_content

    def setup_new_file(self, option_name, basename):
        file = os.path.join(self.kwargs["output_dir"], option_name.replace(":", "_"), str(self.file_counter), basename)
        if self.kwargs["relative_paths"]:
            file = os.path.relpath(file, start=os.curdir)
        self.file_counter += 1
        if not self.kwargs["dry_run"]:
            os.makedirs(os.path.dirname(file), exist_ok=False)

        return file

    def write_file(self, option_name, basename, file_content, mode="w"):
        new_file = self.setup_new_file(option_name, basename)
        if not self.kwargs["dry_run"]:
            with open(new_file, mode) as f:
                f.write(file_content)

        return new_file

    def write_text_file(self, *args):
        return self.write_file(*args, mode="w")

    def write_binary_file(self, *args):
        return self.write_file(*args, mode="wb")

    def get_options_dict(self):
        encoded_options_dict = self.get_encoded_options_dict()
        options_dict = collections.defaultdict(list)
        for option_name, values in encoded_options_dict.items():
            if self.kwargs["skip_corrections"] and option_name.startswith(self._corrections_option_name):
                continue
            for value in values:
                if isinstance(value, dict):
                    for basename, data in value.items():
                        assert isinstance(data, str), self._malformed_init_file_error_msg
                        if data.startswith(PyRosettaInitFileSerializer._prefix_string):
                            file_content = self.format_decode_string(data, option_name)
                            filename = self.write_text_file(option_name, basename, file_content)
                            options_dict[option_name].append(filename)
                        elif data.startswith(PyRosettaInitFileSerializer._prefix_binary):
                            file_content = self.format_decode_binary(data)
                            filename = self.write_binary_file(option_name, basename, file_content)
                            options_dict[option_name].append(filename)
                        else:
                            options_dict[option_name].append(data)
                elif isinstance(value, str):
                    options_dict[option_name].append(value)
                else:
                    raise ValueError(self._malformed_init_file_error_msg)

        return dict(options_dict)

    def get_options(self):
        options_dict = self.get_options_dict()

        return " ".join(
            [
                "-{0} {1}".format(option_name, " ".join(values))
                for option_name, values in options_dict.items()
            ]
        )

    def print_results(self):
        if self.kwargs["dry_run"]:
            print("Dry run PyRosetta initialization from file: {0}".format(self.init_file))
            if self.file_counter == 0:
                print("No PyRosetta input files to decompress.")
            else:
                print("Decompressed {0} PyRosetta input file(s).".format(self.file_counter))
        else:
            print("Initializing PyRosetta from file: {0}".format(self.init_file))
            if self.file_counter == 0:
                print("No PyRosetta input files to decompress.")
            else:
                print("Decompressed {0} PyRosetta input file(s) written to: {1}".format(
                        self.file_counter, self.kwargs["output_dir"]
                    )
                )
        print(
            "Author(s): {0}".format(self.init_dict["author"]),
            "E-mail(s): {0}".format(self.init_dict["email"]),
            "License(s): {0}".format(self.init_dict["license"]),
            "Metadata: {0}".format(self.init_dict["metadata"]),
            "PyRosetta build: {0}".format(self.init_dict["pyrosetta_build"]),
            "Date/Time created (UTC): {0}".format(
                datetime.datetime.strptime(
                    self.init_dict["datetime"],
                    self._strftime_format,
                ).strftime("%b %d, %Y at %I:%M:%S %p")
            ),
            sep=os.linesep,
        )

    def pprint_options(self, options):
        if self.kwargs["dry_run"]:
            print("PyRosetta initialization options from dry run:")
        else:
            print("PyRosetta initialization options:")
        pprint(options)
        if self.kwargs["dry_run"]:
            print("Skipping PyRosetta initialization...")
        else:
            print("Running PyRosetta initialization...")

    def init(self):
        self.validate_init_was_not_called()
        original_pyrosetta_build = self.init_dict["pyrosetta_build"]
        current_pyrosetta_build = self.get_pyrosetta_build()
        if original_pyrosetta_build != current_pyrosetta_build:
            self.pyrosetta_build_warning(original_pyrosetta_build, current_pyrosetta_build)
            if self.kwargs["skip_corrections"] is None:
                self.kwargs["skip_corrections"] = False
        else:
            if self.kwargs["skip_corrections"] is None:
                self.kwargs["skip_corrections"] = True
        options = self.get_options()
        if self.kwargs["verbose"]:
            self.print_results()
            self.pprint_options(options)
        if not self.kwargs["dry_run"]:
            pyrosetta.init(
                options=options,
                extra_options="",
                set_logging_handler=self.kwargs["set_logging_handler"],
                notebook=self.kwargs["notebook"],
                silent=self.kwargs["silent"],
            )


class PyRosettaInitFileParser(object):
    @staticmethod
    def init_from_file(
        init_file,
        dry_run=None,
        output_dir=None,
        skip_corrections=None,
        relative_paths=None,
        max_decompressed_bytes=None,
        database=None,
        verbose=None,
        set_logging_handler=None,
        notebook=None,
        silent=None,
    ):
        """
        Initialize PyRosetta from a '.init' file.

        This method decompresses PyRosetta initialization input files from an input '.init' file into an output directory, and
        then runs `pyrosetta.init` with the cached Rosetta command line options pointing to the files written to the output directory.
        Therefore, it may be helpful to enable the 'dry_run' keyword argument to first inspect the Rosetta command line options
        before committing to writing all PyRosetta input files to disk and running PyRosetta initialization.

        Args:
            init_file: A required `str` object representing the input '.init' file.

        Keyword Args:
            dry_run: An optional `bool` object specifying whether or not to write PyRosetta input files and perform PyRosetta
                initialization. If `True`, then only print the PyRosetta initialization options that would be run if it were `False`.
                Default: False
            output_dir: An optional `str` object representing the output directory into which to decompress PyRosetta input files.
                Default: `./pyrosetta_init_files`
            skip_corrections: An optional `bool` object specifying whether or not to skip any ScoreFunction corrections specified
                in the input '.init' file, which are set in-code upon PyRosetta initialization. If a `NoneType` object is provided,
                then the input ScoreFunction corrections are automatically used for PyRosetta initialization only if the PyRosetta
                version from the '.init' file does not match the current PyRosetta version.
                Default: None
            relative_paths: An optional `bool` object specifying whether or not to initialize PyRosetta with the relative paths
                (with respect to the current working directory) of the files written to the 'output_dir' keyword argument parameter.
                Default: False
            max_decompressed_bytes: An optional `int` object specifying the maximum permitted number of bytes per decompressed PyRosetta
                input file (with a default of 200 MB). If a PyRosetta input file in the input '.init' file exceeds this buffer size
                upon decompression, then a `BufferError` is intentionally raised as a precaution.
                Default: 200_000_000
            database: An optional `str` object representing the path to the Rosetta database. By default, the Rosetta database
                is found using `pyrosetta._rosetta_database_from_env()`, but if the search fails then the Rosetta database path
                may be manually input here.
                Default: None
            verbose: An optional `bool` object specifying whether or not to print PyRosetta initialization information.
                Default: True
            set_logging_handler: An optional object passed to `pyrosetta.init(set_logging_handler=...)` during PyRosetta initialization.
                If `None`, then the default `pyrosetta.init` keyword argument parameter is used. 
                Default: None
            notebook: An optional object passed to `pyrosetta.init(notebook=...)` during PyRosetta initialization.
                If `None`, then the default `pyrosetta.init` keyword argument parameter is used.
                Default: None
            silent: An optional object passed to `pyrosetta.init(silent=...)` during PyRosetta initialization.
                If `None`, then the default `pyrosetta.init` keyword argument parameter is used.
                Default: None

        Returns:
            None
        """
        return PyRosettaInitFileReader(
            init_file,
            dry_run=dry_run,
            output_dir=output_dir,
            skip_corrections=skip_corrections,
            relative_paths=relative_paths,
            max_decompressed_bytes=max_decompressed_bytes,
            database=database,
            verbose=verbose,
            set_logging_handler=set_logging_handler,
            notebook=notebook,
            silent=silent,
        ).init()

    @staticmethod
    def get_init_options_from_file(
        init_file,
        dry_run=True,
        output_dir=None,
        relative_paths=None,
        max_decompressed_bytes=None,
        database=None,
        as_dict=False,
    ):
        """
        Get PyRosetta initialization options from a '.init' file.

        This method returns the PyRosetta initialization options from an input '.init' file without running `pyrosetta.init`. The 
        'dry_run' keyword argument is enabled by default in order to inspect the Rosetta command line options before committing
        to decompressing and writing all PyRosetta input files to disk. If 'dry_run' is disabled, then also decompress the PyRosetta
        input files into an output directory given by the 'output_dir' keyword argument parameter without running `pyrosetta.init`.

        Args:
            init_file: A required `str` object representing the input '.init' file.

        Keyword Args:
            dry_run: An optional `bool` object specifying whether or not to write PyRosetta input files to disk. If `True`, then the
                PyRosetta input files will not be written to disk so the returned options can be inspected (or input manually into 
                `pyrosetta.init` if options do not contain input files).
                Default: True
            output_dir: An optional `str` object representing the output directory into which to decompress PyRosetta input files if
                the 'dry_run' keyword argument parameter is `False`.
                Default: `./pyrosetta_init_files`
            relative_paths: An optional `bool` object specifying whether or not to return the relative paths (with respect to
                the current working directory) of the files written to the 'output_dir' keyword argument parameter.
                Default: False
            max_decompressed_bytes: An optional `int` object specifying the maximum permitted number of bytes per decompressed PyRosetta
                input file (with a default of 200 MB). If a PyRosetta input file in the input '.init' file exceeds this buffer size
                upon decompression, then a `BufferError` is intentionally raised as a precaution.
                Default: 200_000_000
            database: An optional `str` object representing the path to the Rosetta database. By default, the Rosetta database
                is found using `pyrosetta._rosetta_database_from_env()`, but if the search fails then the Rosetta database path
                may be manually input here.
                Default: None
            as_dict: An optional `bool` object specifying whether or not to return the PyRosetta initialization options as a `dict`
                object, otherwise options are returned as a `str` object.
                Default: False

        Raises:
            `ValueError` when the 'as_dict' keyword argument parameter is not a `bool` object.

        Returns:
            A `str` or `dict` object representing the PyRosetta initialization options.
        """
        if not isinstance(as_dict, bool):
            raise ValueError(f"The 'as_dict' keyword argument parameter must be a `bool` object. Received: {type(as_dict)}")
        reader = PyRosettaInitFileReader(
            init_file,
            dry_run=dry_run,
            output_dir=output_dir,
            skip_corrections=False,
            relative_paths=relative_paths,
            max_decompressed_bytes=max_decompressed_bytes,
            database=database,
            verbose=False,
            set_logging_handler=None,
            notebook=None,
            silent=None,
        )
        if as_dict:
            return reader.get_options_dict()
        else:
            return reader.get_options()

    @staticmethod
    def dump_init_file(
        output_filename,
        poses=None,
        author=None,
        email=None,
        license=None,
        metadata=None,
        overwrite=None,
        dry_run=None,
        verbose=None,
    ):
        """
        Write a PyRosetta initialization '.init' file.

        This method uses the `ProtocolSettingsMetric` SimpleMetric to get Rosetta command line options and compresses any PyRosetta input
        files (including subfiles within files) into the output '.init' file. The Rosetta database directory is automatically excluded.
        Only the relative paths of any input directories (from the current working directory) are saved in the Rosetta command line options
        (e.g., '-in:path:bcl /path/to/current/directory/bcl_rosetta' is saved as '-in:path:bcl ./bcl_rosetta'). Therefore, it may be helpful
        to add comments to the 'metadata' keyword argument parameter about specific PyRosetta initialization requirements. PyRosetta
        initialization input files are automatically detected and compressed into the provided 'output_filename' argument parameter,
        and so it can be useful to start with the `dry_run` keyword argument enabled to confirm that the PyRosetta input files are correct.
        Note that automatic detection of PyRosetta input files containing any spaces (e.g., ' ') in file paths or filenames is not supported.

        Args:
            output_filename: A required `str` object representing the output '.init' file.

        Keyword Args:
            poses: An optional `Pose`, `PackedPose`, or iterable of `Pose` or `PackedPose` objects to cache in the output '.init' file.
                Default: []
            author: An optional `str` object representing the author's/authors' name(s) or username(s).
                Default: ""
            email: An optional `str` object representing the author's/authors' email address(es).
                Default: ""
            license: An optional `str` object representing the license(s) for the output '.init' file.
                Default: ""
            metadata: An optional JSON-serializable object representing any additional metadata to save to the output '.init' file.
                Default: {}
            overwrite: An optional `bool` object specifying whether or not to overwrite the output '.init' file if it exists.
                If `False`, then raise an error if the output '.init' file already exists.
                Default: False
            dry_run: An optional `bool` object specifying whether or not to dump the output '.init' file. If `True`, then only print
                the files that would be compressed into the '.init' file if it were `False`.
                Default: False
            verbose: An optional `bool` object specifying whether or not to print PyRosetta '.init' file information.
                Default: True

        Returns:
            None
        """
        return PyRosettaInitFileWriter(
            output_filename,
            poses=poses,
            author=author,
            email=email,
            license=license,
            metadata=metadata,
            overwrite=overwrite,
            dry_run=dry_run,
            verbose=verbose,
        ).dump()

    @staticmethod
    def get_init_options(compressed=False, as_dict=False):
        """
        Get the currently initialized Rosetta command line options using the `ProtocolSettingsMetric` SimpleMetric, including the
        Rosetta database.

        Keyword Args:
            compressed: An optional `bool` object specifying whether or not to compress any input files (including files containing
                lists of files) in memory, and return only the relative paths of any input directories (from the current working
                directory) in the Rosetta command line options (e.g., '-in:path:bcl /path/to/current/directory/bcl_rosetta' is
                returned as '-in:path:bcl ./bcl_rosetta').
                Default: False
            as_dict: An optional `bool` object specifying whether or not to return the PyRosetta initialization options as a `dict`
                object, otherwise options are returned as a `str` object.
                Default: False

        Raises:
            `ValueError` when the 'compressed' or 'as_dict' keyword argument parameters are not `bool` objects.
            `NotImplementedError` when `compressed=True` and `as_dict=False`.

        Returns:
            A `str` or `dict` object representing the PyRosetta initialization options.
        """
        if not isinstance(compressed, bool):
            raise ValueError(f"The 'compressed' keyword argument parameter must be a `bool` object. Received: {type(compressed)}")
        if not isinstance(as_dict, bool):
            raise ValueError(f"The 'as_dict' keyword argument parameter must be a `bool` object. Received: {type(as_dict)}")
        with tempfile.TemporaryDirectory() as tmp_dir:
            writer = PyRosettaInitFileWriter(
                os.path.join(tmp_dir, "tmp.init"),
                poses=None,
                author=None,
                email=None,
                license=None,
                metadata=None,
                overwrite=False,
                dry_run=False,
                verbose=False,
            )
        if compressed:
            if as_dict:
                return writer.get_encoded_options_dict()
            else:
                raise NotImplementedError("Formatting compressed PyRosetta initialization options into a `str` object is not supported.")
        else:
            if as_dict:
                return writer.get_options_dict()
            else:
                return writer.get_options_str()
