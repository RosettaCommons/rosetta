# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import pandas
    import toolz
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.io' requires the "
        + "third-party packages 'pandas' and 'toolz' as a dependencies!\n"
        + "Please install these packages into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/toolz/\n"
        + "https://pypi.org/project/pandas/\n"
    )
    raise

import bz2
import collections
import json
import logging
import os
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import re
import uuid
import warnings

from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime
from pyrosetta.rosetta.core.pose import Pose, add_comment, get_all_comments
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.exceptions import PyRosettaIsNotInitializedError
from pyrosetta.rosetta.basic import was_init_called
from pyrosetta.secure_unpickle import SecureSerializerBase
from pyrosetta.utility.initialization import (
    PyRosettaInitDictWriter,
    PyRosettaInitFileReader,
)
from typing import (
    Any,
    Dict,
    Generic,
    Iterable,
    List,
    NoReturn,
    Optional,
    Tuple,
    TypeVar,
    Union,
)
from urllib.parse import urlparse, urlunparse

from pyrosetta.distributed.cluster.config import source_domains
from pyrosetta.distributed.cluster.exceptions import OutputError
from pyrosetta.distributed.cluster.init_files import InitFileSigner
from pyrosetta.distributed.cluster.logging_support import RedirectToLogger
from pyrosetta.distributed.cluster.serialization import Serialization, update_scores


METADATA_INPUT_DECOY_KEY: str = "idx_poses_input"
METADATA_OUTPUT_DECOY_KEY: str = "idx_poses_output"

G = TypeVar("G")


class IO(Generic[G]):
    """Input/Output methods for PyRosettaCluster."""

    DATETIME_FORMAT: str = "%Y-%m-%d %H:%M:%S.%f"
    REMARK_FORMAT: str = "REMARK PyRosettaCluster: "

    def _get_instance_and_metadata(
        self, kwargs: Dict[Any, Any]
    ) -> Tuple[Dict[Any, Any], Dict[Any, Any]]:
        """
        Get the current state of the PyRosettaCluster instance, and split the
        kwargs into the PyRosettaCluster instance kwargs and ancillary metadata.
        """

        instance_state = self.__getstate__()
        if isinstance(instance_state, tuple):
            instance_state = dict(zip(self.__slots__, self.__getstate__()))
        assert isinstance(instance_state, dict)
        instance_state.pop("security", None)
        instance_state.pop("client", None)
        instance_state.pop("clients", None)
        instance_state.pop("output_init_file", None)
        instance_state.pop("logging_address", None)
        instance_kwargs = self.serializer.deepcopy_kwargs(instance_state)
        for i in self.__attrs_attrs__:
            if not i.init:
                instance_kwargs.pop(i.name)
            if i.name == "input_packed_pose":
                instance_kwargs.pop(i.name, None)
        instance_kwargs["tasks"] = kwargs.pop("task")
        instance_kwargs["seeds"] = kwargs.pop("seeds")
        instance_kwargs["decoy_ids"] = kwargs.pop("decoy_ids")

        return instance_kwargs, kwargs

    def _get_output_dir(self, decoy_dir: str) -> str:
        """Get the output directory in which to write files to disk."""

        zfill_value = 4
        max_dir_depth = 1000
        decoy_dir_list = os.listdir(decoy_dir)
        if not decoy_dir_list:
            new_dir = str(0).zfill(zfill_value)
            output_dir = os.path.join(decoy_dir, new_dir)
            os.mkdir(output_dir)
        else:
            top_dir = list(reversed(sorted(decoy_dir_list)))[0]
            if len(os.listdir(os.path.join(decoy_dir, top_dir))) < max_dir_depth:
                output_dir = os.path.join(decoy_dir, top_dir)
            else:
                new_dir = str(int(top_dir) + 1).zfill(zfill_value)
                output_dir = os.path.join(decoy_dir, new_dir)
                os.mkdir(output_dir)

        return output_dir

    @staticmethod
    def _filter_scores_dict(scores_dict: Dict[Any, Any]) -> Dict[Any, Any]:
        for key in list(scores_dict.keys()):
            try:
                IO._dump_json(scores_dict[key])
            except:
                logging.warning(
                    f"Removing score key '{key}' with value of type '{type(scores_dict[key])}' before "
                    + "saving PyRosettaCluster result! Only JSON-serializable score values can be written to "
                    + "output files. Consider custom serializing the value to save this score or removing the "
                    + "key from the `pose.cache` dictionary to remove this warning message."
                )
                scores_dict.pop(key, None)

        return scores_dict

    def _format_result(self, result: Union[Pose, PackedPose]) -> Tuple[str, Dict[Any, Any], PackedPose]:
        """
        Given a `Pose` or `PackedPose` object, return a tuple containing
        the pdb string and a scores dictionary.
        """

        _pdbstring = io.to_pdbstring(result)
        _scores_dict = update_scores(PackedPose(result)).scores
        _filtered_scores_dict = IO._filter_scores_dict(self.serializer.deepcopy_kwargs(_scores_dict))

        return (result, _pdbstring, _scores_dict, _filtered_scores_dict)

    def _parse_results(
        self,
        results: Union[
            Iterable[Optional[Union[Pose, PackedPose, bytes]]],
            Optional[Union[Pose, PackedPose]],
        ],
    ) -> Union[List[Tuple[str, Dict[Any, Any]]], NoReturn]:
        """
        Format output results on distributed worker. Input argument `results` can be a
        `Pose`, `PackedPose`, or `None` object, or a `list` or `tuple` of `Pose` and/or `PackedPose`
        objects, or an empty `list` or `tuple`. Returns a list of tuples, each tuple
        containing the pdb string and a scores dictionary.
        """
        if isinstance(results, bytes):
            results = self.serializer.decompress_packed_pose(results)

        if isinstance(results, (Pose, PackedPose)):
            if not io.to_pose(results).empty():
                out = [self._format_result(results)]
            else:
                out = []
        elif isinstance(results, collections.abc.Iterable):
            out = []
            for result in results:
                if isinstance(results, bytes):
                    result = self.serializer.decompress_packed_pose(result)
                if isinstance(result, (Pose, PackedPose)):
                    if not io.to_pose(result).empty():
                        out.append(self._format_result(result))
                else:
                    raise OutputError(result)
        elif not results:
            out = []
        else:
            raise OutputError(results)

        return out

    def _process_kwargs(self, kwargs: Dict[Any, Any]) -> Dict[Any, Any]:
        """
        Remove seed specification from 'extra_options' or 'options',
        and remove protocols_key from kwargs.
        """

        for k in list(kwargs.keys()):
            if not k.startswith("PyRosettaCluster_"):
                kwargs.pop(k, None)
        kwargs.pop(self.protocols_key, None)
        kwargs.pop("PyRosettaCluster_protocol_number", None)
        kwargs.pop("PyRosettaCluster_protocol_name", None)
        kwargs.pop("PyRosettaCluster_seed", None)
        kwargs.pop("PyRosettaCluster_output_path", None)
        kwargs.pop("PyRosettaCluster_tmp_path", None)
        kwargs.pop("PyRosettaCluster_client_repr", None)
        task_protocols = kwargs["PyRosettaCluster_protocols"]
        task_seeds = kwargs["PyRosettaCluster_seeds"]
        task_decoy_ids = kwargs["PyRosettaCluster_decoy_ids"]
        assert len(task_protocols) == len(
            task_seeds
        ), "Task protocols and seeds must be the same length."
        assert len(task_protocols) == len(
            task_decoy_ids
        ), "Task protocols and decoy_ids must be the same length."
        for i in range(len(task_protocols)):
            assert (
                task_protocols[i] == task_seeds[i][0]
            ), "Task protocol elements and seed keys must be the same."
            assert (
                task_protocols[i] == task_decoy_ids[i][0]
            ), "Task protocol elements and decoy_id keys must be the same."
        kwargs["PyRosettaCluster_seeds"] = [v for (k, v) in task_seeds]
        kwargs["PyRosettaCluster_decoy_ids"] = [v for (k, v) in task_decoy_ids]

        return kwargs

    def _get_init_file_json(self, packed_pose: PackedPose) -> str:
        """Return a PyRosetta initialization file as a JSON-serialized string."""

        metadata, poses = sign_init_file_metadata_and_poses(
            input_packed_pose=self.input_packed_pose,
            output_packed_pose=packed_pose,
        )
        writer = PyRosettaInitDictWriter(
            poses=poses,
            author=self.author,
            email=self.email,
            license=self.license,
            metadata=metadata,
            overwrite=False,
            dry_run=self.dry_run,
            verbose=False,
        )

        return writer.get_json()

    @staticmethod
    def _add_pose_comment(packed_pose: PackedPose, pdbfile_data: str) -> PackedPose:
        """Cache simulation data as a pose comment."""

        _pose = packed_pose.pose.clone()
        add_comment(
            _pose,
            IO.REMARK_FORMAT.rstrip(), # Remove extra space since `add_comment` adds a space
            pdbfile_data,
        )

        return io.to_packed(_pose)

    @staticmethod
    def _dump_json(data: Dict[str, Any]) -> str:
        """Return JSON-serialized data."""

        return json.dumps(
            data,
            skipkeys=False,
            ensure_ascii=False,
            check_circular=True,
            allow_nan=False,
            cls=None,
            indent=None,
            separators=(", ", ": "),
            default=None,
            sort_keys=False,
        )

    def _save_results(self, results: Any, kwargs: Dict[Any, Any]) -> None:
        """Write results and kwargs to disk."""

        if self.dry_run:
            logging.info(
                "PyRosettaCluster `dry_run` attribute is enabled. Skipping saving!"
            )
            return

        # Parse and save results
        for packed_pose, pdbstring, scores, scores_json in self._parse_results(results):
            kwargs = self._process_kwargs(kwargs)
            output_dir = self._get_output_dir(decoy_dir=self.decoy_path)
            decoy_name = "_".join([self.simulation_name, uuid.uuid4().hex])
            output_file_extension = self.output_decoy_types[0]
            output_file = os.path.join(output_dir, decoy_name + output_file_extension)
            if self.compressed:
                output_file += ".bz2"
            extra_kwargs = {
                "PyRosettaCluster_decoy_name": decoy_name,
                "PyRosettaCluster_output_file": output_file,
            }
            extra_kwargs["PyRosettaCluster_environment_manager"] = self.environment_manager
            if self.manifest:
                extra_kwargs["PyRosettaCluster_manifest"] = self.manifest
            if self.manifest_filename:
                extra_kwargs["PyRosettaCluster_manifest_filename"] = self.manifest_filename
            if os.path.isfile(self.environment_file):
                extra_kwargs["PyRosettaCluster_environment_file"] = self.environment_file
            if os.path.isfile(self.output_init_file):
                extra_kwargs["PyRosettaCluster_init_file"] = self.output_init_file
            if "PyRosettaCluster_datetime_start" in kwargs:
                datetime_end = datetime.now().strftime(self.DATETIME_FORMAT)
                duration = str(
                    (
                        datetime.strptime(datetime_end, self.DATETIME_FORMAT)
                        - datetime.strptime(
                            kwargs["PyRosettaCluster_datetime_start"],
                            self.DATETIME_FORMAT,
                        )
                    ).total_seconds()
                )  # For build-in functions
                extra_kwargs.update(
                    {
                        "PyRosettaCluster_datetime_end": datetime_end,
                        "PyRosettaCluster_total_seconds": duration,
                    }
                )
            instance, metadata = self._get_instance_and_metadata(
                toolz.dicttoolz.keymap(
                    lambda k: k.split("PyRosettaCluster_")[-1],
                    toolz.dicttoolz.merge(extra_kwargs, kwargs),
                )
            )
            instance_metadata = {
                "instance": collections.OrderedDict(sorted(instance.items())),
                "metadata": collections.OrderedDict(sorted(metadata.items())),
            }
            simulation_data = self.serializer.deepcopy_kwargs(instance_metadata)
            simulation_data["scores"] = collections.OrderedDict(sorted(scores.items()))
            simulation_data_json = self.serializer.deepcopy_kwargs(instance_metadata)
            simulation_data_json["scores"] = collections.OrderedDict(sorted(scores_json.items()))
            pdbfile_data = IO._dump_json(simulation_data_json)
            # Output PDB file
            if ".pdb" in self.output_decoy_types:
                # Write full .pdb record
                pdbstring_data = pdbstring + os.linesep + self.REMARK_FORMAT + pdbfile_data
                output_pdb_file = os.path.join(output_dir, decoy_name + ".pdb")
                if self.compressed:
                    output_pdb_file += ".bz2"
                    with open(output_pdb_file, "wb") as f:
                        f.write(bz2.compress(str.encode(pdbstring_data)))
                else:
                    with open(output_pdb_file, "w") as f:
                        f.write(pdbstring_data)

            # Output pickled Pose file
            if ".pkl_pose" in self.output_decoy_types:
                _packed_pose = IO._add_pose_comment(packed_pose, pdbfile_data)
                output_pkl_pose_file = os.path.join(output_dir, decoy_name + ".pkl_pose")
                if self.compressed:
                    output_pkl_pose_file += ".bz2"
                    with open(output_pkl_pose_file, "wb") as f:
                        f.write(bz2.compress(io.to_pickle(_packed_pose)))
                else:
                    with open(output_pkl_pose_file, "wb") as f:
                        f.write(io.to_pickle(_packed_pose))

            # Output base64-encoded pickled Pose file
            if ".b64_pose" in self.output_decoy_types:
                _packed_pose = IO._add_pose_comment(packed_pose, pdbfile_data)
                output_b64_pose_file = os.path.join(output_dir, decoy_name + ".b64_pose")
                if self.compressed:
                    output_b64_pose_file += ".bz2"
                    with open(output_b64_pose_file, "wb") as f:
                        f.write(bz2.compress(str.encode(io.to_base64(_packed_pose))))
                else:
                    with open(output_b64_pose_file, "w") as f:
                        f.write(io.to_base64(_packed_pose))

            # Output PyRosetta initialization file
            if ".init" in self.output_decoy_types:
                _packed_pose = IO._add_pose_comment(packed_pose, pdbfile_data)
                init_file_json = self._get_init_file_json(_packed_pose)
                output_init_file = os.path.join(output_dir, decoy_name + ".init")
                if self.compressed:
                    output_init_file += ".bz2"
                    with open(output_init_file, "wb") as f:
                        f.write(bz2.compress(str.encode(init_file_json)))
                else:
                    with open(output_init_file, "w") as f:
                        f.write(init_file_json)

            # Output JSON-encoded scorefile
            if ".json" in self.output_scorefile_types:
                if self.simulation_records_in_scorefile:
                    scorefile_data = pdbfile_data
                else:
                    scorefile_data = IO._dump_json(
                        {
                            metadata["output_file"]: collections.OrderedDict(
                                sorted(scores_json.items())
                            ),
                        }
                    )
                # Append data to scorefile
                with open(self.scorefile_path, "a") as f:
                    f.write(scorefile_data + os.linesep)

            # Output pickled `pandas.DataFrame` scorefile
            for extension in self.output_scorefile_types:
                if extension != ".json":
                    _scorefile_path = os.path.splitext(self.scorefile_path)[0] + extension
                    if self.simulation_records_in_scorefile:
                        _scorefile_data = {
                            metadata["output_file"]: collections.OrderedDict(simulation_data),
                        }
                    else:
                        _scorefile_data = {
                            metadata["output_file"]: simulation_data["scores"]
                        }
                    df = pandas.DataFrame().from_dict(_scorefile_data, orient="index")
                    # Append data to scorefile
                    if os.path.isfile(_scorefile_path):
                        df_chunk = secure_read_pickle(_scorefile_path, compression="infer")
                        df = pandas.concat([df_chunk, df])
                    df.to_pickle(_scorefile_path, compression="infer", protocol=SecureSerializerBase._pickle_protocol)

    def _cache_manifest(self) -> None:
        """Cache the pixi manifest 'pixi.toml' or 'pyproject.toml' file string and manifest format."""

        if self.environment_manager == "pixi":
            # https://pixi.sh/dev/reference/environment_variables/#environment-variables-set-by-pixi
            manifest_path = os.environ.get("PIXI_PROJECT_MANIFEST", "")
            if manifest_path:
                with open(manifest_path, "r") as f:
                    self.manifest = sanitize_urls(f.read())
            else:
                # https://pixi.sh/dev/python/tutorial/#pixitoml-and-pyprojecttoml
                for filename in ("pixi.toml", "pyproject.toml"):
                    manifest_path = os.path.join(os.getcwd(), filename)
                    if os.path.isfile(manifest_path):
                        with open(manifest_path, "r") as f:
                            self.manifest = sanitize_urls(f.read())
                        break
                else:
                    self.manifest = ""

            # Cache TOML filename as manifest format
            if self.manifest:
                self.manifest_filename = os.path.basename(manifest_path)
            else:
                logging.warning(
                    (
                        "PyRosettaCluster could not detect the pixi manifest file! "
                        "It is recommended to commit the pixi manifest file to the "
                        "git repository to reproduce the pixi project later."
                    ),
                    UserWarning,
                    stacklevel=3,
                )
        else:
            self.manifest = ""
            self.manifest_filename = ""

    def _write_environment_file(self, filename: str) -> None:
        """
        Write the conda/mamba YML, uv requirements, or pixi lock file string to the input filename.
        If pixi is used as the environment manager, also write the manifest file string to a separate filename.
        """

        if (
            (not self.simulation_records_in_scorefile)
            and (not self.dry_run)
            and self.environment
        ):
            with open(filename, "w") as f:
                f.write(self.environment)

            if self.environment_manager == "pixi" and self.manifest and self.manifest_filename:
                toml_file = "_".join(filename.split("_")[:-1] + [self.manifest_filename])
                with open(toml_file, "w") as f:
                    f.write(self.manifest)

    def _write_init_file(self) -> None:
        """Maybe write PyRosetta initialization file to the input filename."""

        if self.output_init_file != "":
            if self.compressed and self.output_init_file.endswith(".init"):
                self.output_init_file += ".bz2"
            if not self.output_init_file.endswith((".init", ".init.bz2")):
                raise ValueError(
                    "The output PyRosetta initialization file must end in '.init' or '.init.bz2'. "
                    + "The current PyRosettaCluster 'output_init_file' instance attribute is: "
                    + f"'{self.output_init_file}'"
                )
            self._dump_init_file(
                self.output_init_file,
                input_packed_pose=self.input_packed_pose,
                output_packed_pose=None,
                verbose=True,
            )

    def _dump_init_file(
        self,
        filename: str,
        input_packed_pose: Optional[PackedPose] = None,
        output_packed_pose: Optional[PackedPose] = None,
        verbose: bool = True,
    ) -> None:
        """Dump compressed PyRosetta initialization input files and poses to the input filename."""

        out = RedirectToLogger(logging.INFO)
        err = RedirectToLogger(logging.ERROR)
        with redirect_stdout(out), redirect_stderr(err):
            metadata, poses = sign_init_file_metadata_and_poses(
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
            )
            dump_init_file_kwargs = dict(
                poses=poses,
                author=self.author,
                email=self.email,
                license=self.license,
                metadata=metadata,
                overwrite=True,
                dry_run=self.dry_run,
                verbose=verbose,
            )
            _logging_debug_msg = "Successfully ran `{0}` on the host node."
            _logging_error_msg = (
                "{0}: {1}. `{2}` did not run successfully, so PyRosetta initialization input data may not "
                "be saved! It is recommended to run `pyrosetta.dump_init_file()` to reproduce PyRosetta "
                "initialization options on the host node later."
            )
            if self.compressed:
                try:
                    writer = PyRosettaInitDictWriter(**dump_init_file_kwargs)
                    init_file_json = writer.get_json()
                    if verbose:
                        writer.print_cached_files(filename, dump_init_file_kwargs["dry_run"])
                    with open(filename, "wb") as f:
                        f.write(bz2.compress(str.encode(init_file_json)))
                    logging.debug(_logging_debug_msg.format("PyRosettaInitDictWriter().get_json()"))
                except Exception as ex:
                    logging.error(_logging_error_msg.format(type(ex).__name__, ex, "PyRosettaInitDictWriter().get_json()"))
            else:
                try:
                    pyrosetta.dump_init_file(filename, **dump_init_file_kwargs)
                    logging.debug(_logging_debug_msg.format("pyrosetta.dump_init_file()"))
                except Exception as ex:
                    logging.error(_logging_error_msg.format(type(ex).__name__, ex, "pyrosetta.dump_init_file()"))
        out.flush()
        err.flush()


def verify_init_file(
    init_file: str,
    input_packed_pose: Optional[PackedPose],
    output_packed_pose: Optional[PackedPose],
    metadata: Dict[str, Any],
) -> Optional[NoReturn]:
    """Verify a PyRosetta initialization file."""

    @toolz.functoolz.curry
    def _verify_signer(_signer: InitFileSigner, _sha256: str, _signature: str) -> Optional[NoReturn]:
        """Verify that current PyRosetta and PyRosettaCluster versions match that dumped in the '.init' file"""

        _init_file_err_msg = (
            f"Failed to verify data integrity of the input PyRosetta initialization file: '{init_file}'"
        )
        _err_msg = (
            "The expected {0} differs from the metadata '{1}' key value in the '.init' file! {2}. The "
            + "simulation cannot necessarily be reproduced! To override '.init' file data verification, "
            + "delete the '{1}' key and value from the metadata in the '.init' file.\n"
            + "Expected: '{3}'\n"
            + "Value:    '{4}'\n"
            + _init_file_err_msg
        )
        _sha256_err_msg = (
            "Either (1) the current PyRosetta build differs from that used to write the original '.init' file, "
            + f"or (2) the '.init' file 'poses' key value (or the metadata key '{METADATA_INPUT_DECOY_KEY}' and "
            + f"'{METADATA_OUTPUT_DECOY_KEY}' and their values) has been altered from the original simulation or "
            + "'.init' file export"
        )
        _signature_err_msg = (
            "Either (1) the current PyRosetta build differs from that used to write the original '.init' file, (2) "
            + "the '.init' file metadata does not match the expected output format by PyRosettaCluster, or (3) the "
            + "'.init' file data has been altered from the original PyRosettaCluster simulation or '.init' file export"
        )

        if (_sha256 is not None) and (not _signer.verify_sha256(_sha256)):
            raise ValueError(
                _err_msg.format(
                    "SHA256 value",
                    "sha256",
                    _sha256_err_msg,
                    _signer.sign_sha256(),
                    _sha256,
                )
            )
        if (_signature is not None) and (not _signer.verify_signature(_signature)):
            raise ValueError(
                _err_msg.format(
                    "PyRosettaCluster signature",
                    "signature",
                    _signature_err_msg,
                    _signer.sign_digest(),
                    _signature,
                )
            )
        if all(val is not None for val in (_sha256, _signature)):
            assert _signer.verify(_sha256, _signature), _err_msg


    metadata = Serialization(compression=False).deepcopy_kwargs(metadata)
    verifier = _verify_signer(
        _sha256=metadata.pop("sha256", None),
        _signature=metadata.pop("signature", None),
    )
    _io_kwargs = dict(
        input_packed_pose=input_packed_pose,
        output_packed_pose=output_packed_pose,
    )

    # Verify metadata in the '.init' file matches originally dumped metadata
    verifier(
        _signer=InitFileSigner(
            **toolz.dicttoolz.merge(_io_kwargs, dict(metadata=metadata))
        ),
    )

    # Verify metadata dumped in the '.init' file matches expected metadata format from PyRosettaCluster
    expected_metadata, _poses = sign_init_file_metadata_and_poses(**_io_kwargs)
    expected_metadata.pop("sha256", None)
    expected_metadata.pop("signature", None)
    verifier(
        _signer=InitFileSigner(
            **toolz.dicttoolz.merge(_io_kwargs, dict(metadata=expected_metadata))
        ),
    )


def sign_init_file_metadata_and_poses(
    input_packed_pose: Optional[PackedPose] = None,
    output_packed_pose: Optional[PackedPose] = None,
) -> Tuple[Dict[str, Any], List[PackedPose]]:
    """Sign PyRosetta initialization file 'metadata' and 'poses' keys."""

    metadata = {}
    metadata["comment"] = "Generated by PyRosettaCluster"
    metadata["version"] = pyrosetta.distributed.cluster.__version__
    poses = []
    if isinstance(output_packed_pose, PackedPose): # Set output decoy to front of poses list if present
        _comments = dict(get_all_comments(output_packed_pose.pose))
        _key = IO.REMARK_FORMAT.rstrip() # Remove extra space since `add_comment` adds a space
        if _key not in _comments.keys():
            raise NotImplementedError(
                "Signing '.init' file metadata and poses when PyRosettaCluster simulation data is not "
                + f"cached in the output `PackedPose` object comments is not supported: {_comments}"
            )
        poses.append(output_packed_pose)
        metadata[METADATA_OUTPUT_DECOY_KEY] = poses.index(output_packed_pose)
    if isinstance(input_packed_pose, PackedPose):
        poses.append(input_packed_pose)
        metadata[METADATA_INPUT_DECOY_KEY] = poses.index(input_packed_pose)
    signer = InitFileSigner(
        input_packed_pose=input_packed_pose,
        output_packed_pose=output_packed_pose,
        metadata=metadata,
    )
    metadata.update(signer.sign())

    return metadata, poses


def get_poses_from_init_file(
    init_file: str,
    verify: bool = False,
) -> Union[Tuple[Optional[PackedPose], Optional[PackedPose]], NoReturn]:
    """
    Return a `tuple` object of the input `PackedPose` object and the output `PackedPose` object
    from a '.init' file, and optionally verify PyRosettaCluster metadata in the '.init' file.
    """

    @toolz.functoolz.curry
    def _maybe_to_packed(
        key: str, poses: List[str], metadata: Dict[str, str]
    ) -> Union[Optional[PackedPose], NoReturn]:
        assert isinstance(metadata, dict), (
            "The PyRosetta initialization file 'metadata' key value must be a `dict` object. "
            + f"Received: {type(metadata)}"
        )
        assert isinstance(poses, list), (
            "The PyRosetta initialization file 'poses' key value must be a `list` object. "
            + f"Received: {type(poses)}"
        )
        if key in metadata:
            idx = metadata[key]
            if isinstance(idx, int):
                return io.to_packed(io.to_pose(poses[idx]))
            else:
                raise TypeError(
                    f"The PyRosetta initialization file metadata '{key}' key value must be an `int` object. "
                    + f"Received: {type(idx)}"
                )
        else:
            return None

    if not was_init_called():
        raise PyRosettaIsNotInitializedError(
            f"Please first initialize PyRosetta with the 'init_file' argument parameter: '{init_file}'"
        )
    if not isinstance(init_file, str):
        raise TypeError(f"The 'init_file' argument parameter must be a `str` object. Received: {type(init_file)}")
    if init_file.endswith(".init.bz2"):
        with open(init_file, "rb") as fbz2:
            init_dict = PyRosettaInitFileReader.from_json(bz2.decompress(fbz2.read()).decode())
    elif init_file.endswith(".init"):
        init_dict = PyRosettaInitFileReader.read_json(init_file)
    else:
        raise ValueError("The 'init_file' argument parameter must end with '.init' or '.init.bz2'.")

    metadata = init_dict["metadata"]
    poses = init_dict["poses"]
    _get_packed_pose = _maybe_to_packed(poses=poses, metadata=metadata)
    input_packed_pose = _get_packed_pose(METADATA_INPUT_DECOY_KEY)
    output_packed_pose = _get_packed_pose(METADATA_OUTPUT_DECOY_KEY)

    if verify:
        verify_init_file(init_file, input_packed_pose, output_packed_pose, metadata)

    return (input_packed_pose, output_packed_pose)


def secure_read_pickle(
    filepath_or_buffer: str,
    compression: Optional[Union[str, Dict[str, Any]]] = "infer",
    storage_options: Optional[Dict[str, Any]] = None,
) -> pandas.DataFrame:
    """
    Secure replacement for `pandas.read_pickle()` for file-like objects.

    Usage requires adding 'pandas' as a secure package to unpickle in PyRosetta:
    `pyrosetta.secure_unpickle.add_secure_package('pandas')`
    """
    with pandas.io.common.get_handle(
        filepath_or_buffer,
        "rb",
        compression=compression,
        is_text=False,
        storage_options=storage_options,
    ) as handles:
        return SecureSerializerBase.secure_load(handles.handle)


def sanitize_urls(yml_str: str) -> str:
    """
    Scan the input string and sanitize any URLs that include
    credentials for source domains, returning the updated string.
    """

    def sanitize_url(url: str) -> str:
        """Remove username and password from URLs pointing to source domains."""
        parsed = urlparse(url)

        # No credentials present
        if "@" not in parsed.netloc:
            return url

        # Split credentials from host
        _credentials, host = parsed.netloc.split("@", 1)
        host_domain = host.split(":", 1)[0]  # Remove port if present

        # Only sanitize if the domain is a source domain
        if host_domain not in source_domains:
            return url

        # Build sanitized URL
        sanitized = parsed._replace(netloc=host)
        sanitized_url = urlunparse(sanitized)

        # Warn without leaking credentials
        warnings.warn(
            (
                "PyRosettaCluster automatically removed embedded credentials from the "
                f"conda channel '{host_domain}' while processing the environment file. "
                "These credentials are no longer required by this conda channel. "
                "Please remove them from your configuration to silence this warning."
            ),
            UserWarning,
            stacklevel=2,
        )

        return sanitized_url

    # Match all URLs (i.e., `http://` and `https://` with or without credentials)
    url_regex = re.compile(r'https?://[^\s\'"]+')

    def replacer(match: re.Match) -> str:
        url = match.group(0)
        return sanitize_url(url)

    yml_sanitized_str = url_regex.sub(replacer, yml_str)

    return yml_sanitized_str
