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
        + "third-party packages 'pandas' and 'toolz' as dependencies!\n"
        + "Please install these packages into your virtual environment. "
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

from contextlib import (
    redirect_stdout,
    redirect_stderr,
)
from datetime import datetime
from pyrosetta.rosetta.core.pose import (
    Pose,
    add_comment,
    get_all_comments,
)
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.exceptions import PyRosettaIsNotInitializedError
from pyrosetta.rosetta.basic import was_init_called
from pyrosetta.secure_unpickle import SecureSerializerBase
from pyrosetta.utility.initialization import (
    PyRosettaInitDictWriter,
    PyRosettaInitFileReader,
)
from urllib.parse import (
    urlparse,
    urlunparse,
)

from pyrosetta.distributed.cluster.config import source_domains
from pyrosetta.distributed.cluster.exceptions import OutputError
from pyrosetta.distributed.cluster.init_files import InitFileSigner
from pyrosetta.distributed.cluster.logging_support import RedirectToLogger
from pyrosetta.distributed.cluster.serialization import (
    Serialization,
    update_scores,
)
from pyrosetta.distributed.cluster.type_defs import (
    Any,
    Dict,
    Iterable,
    List,
    Optional,
    PoseOrPackedPose,
    Tuple,
    Union,
)

METADATA_INPUT_DECOY_KEY: str = "idx_poses_input"
METADATA_OUTPUT_DECOY_KEY: str = "idx_poses_output"


class IO:
    """Input/Output methods for `PyRosettaCluster`."""

    DATETIME_FORMAT: str = "%Y-%m-%d %H:%M:%S.%f"
    REMARK_FORMAT: str = "REMARK PyRosettaCluster: "

    def _get_instance_and_metadata(
        self, kwargs: Dict[str, Any]
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """
        Get the current state of the `PyRosettaCluster` instance, and split the input keyword arguments into
        the `PyRosettaCluster` instance attributes and ancillary metadata.
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
        instance_state.pop("max_task_replicas", None)
        instance_state.pop("task_registry", None)
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
    def _filter_scores_dict(scores_dict: Dict[str, Any]) -> Dict[str, Any]:
        """Filter for JSON-serializable scoring data."""

        for key in list(scores_dict.keys()):
            try:
                IO._dump_json(scores_dict[key])
            except:
                logging.warning(
                    f"Removing score key '{key}' with value of type '{type(scores_dict[key])}' before saving "
                    + "`PyRosettaCluster` result! Only JSON-serializable scoring data can be written to output "
                    + "decoy files and scorefiles. Consider custom serializing the value to save this score "
                    + "or removing the key from the `Pose.cache` dictionary to silence this warning message."
                )
                scores_dict.pop(key, None)

        return scores_dict

    def _format_result(
        self, result: PoseOrPackedPose
    ) -> Tuple[PackedPose, str, Dict[str, Any], Dict[str, Any]]:
        """
        Given a `Pose` or `PackedPose` object, return a `tuple` object containing the `Pose` or `PackedPose`
        object, and its PDB string, `Pose.cache` dictionary, and JSON-serializable `Pose.cache` dictionary.

        *Warning*: This method uses the `pickle` module to deserialize pickled `Pose` objects and arbitrary
        Python types in `Pose.cache` dictionary. Using the `pickle` module is not secure, so please only run
        with input files you trust. Learn more about the `pickle` module and its security
        `here <https://docs.python.org/3/library/pickle.html>`_.
        """

        _pdbstring = io.to_pdbstring(result)
        _scores_dict = dict(update_scores(PackedPose(result)).pose.cache)
        _filtered_scores_dict = IO._filter_scores_dict(self.serializer.deepcopy_kwargs(_scores_dict))

        return (result, _pdbstring, _scores_dict, _filtered_scores_dict)

    def _parse_results(
        self,
        results: Optional[Union[bytes, PoseOrPackedPose, Iterable[Union[bytes, PoseOrPackedPose]]]],
    ) -> List[Tuple[str, Dict[str, Any]]]:
        """
        Format output results from a Dask worker.

        *Warning*: This method uses the `pickle` module to deserialize pickled `Pose` objects and arbitrary
        Python types in `Pose.cache` dictionary. Using the `pickle` module is not secure, so please only run
        with input files you trust. Learn more about the `pickle` module and its security
        `here <https://docs.python.org/3/library/pickle.html>`_.

        Args:
            `results`: `Pose | PackedPose | bytes | Iterable[Pose | PackedPose | bytes] | None`
                An `Pose`, `PackedPose`, `bytes` or `None` object, or an iterable of `Pose`, `PackedPose`,
                or `bytes` objects.

        Returns:
            A `list` object of `tuple` objects, where each `tuple` object contains a PDB string, `Pose.cache`
            dictionary, and JSON-serializable `Pose.cache` dictionary.
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

    def _process_kwargs(self, kwargs: Dict[str, Any]) -> Dict[str, Any]:
        """Parse a returned task dictionary."""

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
        """
        Return a PyRosetta initialization file as a JSON-serialized string.

        *Warning*: This method uses the `pickle` module to deserialize pickled `Pose` objects. Using the
        `pickle` module is not secure, so please only run with input files you trust. Learn more about the
        `pickle` module and its security `here <https://docs.python.org/3/library/pickle.html>`_.
        """

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
        """
        Cache simulation data as a `Pose` comment.

        *Warning*: This method uses the `pickle` module to deserialize pickled `Pose` objects. Using the
        `pickle` module is not secure, so please only run with input files you trust. Learn more about the
        `pickle` module and its security `here <https://docs.python.org/3/library/pickle.html>`_.
        """

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

    def _save_results(self, results: Optional[bytes], kwargs: Dict[str, Any]) -> None:
        """
        Write output results to disk.

        *Warning*: This method uses the `pickle` module to deserialize pickled `Pose` objects and arbitrary
        Python types in `Pose.cache` dictionary. Using the `pickle` module is not secure, so please only run
        with input files you trust. Learn more about the `pickle` module and its security
        `here <https://docs.python.org/3/library/pickle.html>`_.
        """

        if self.dry_run:
            logging.info(
                "The `dry_run` instance attribute is set to `True`. Skipping saving!"
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
            if self.toml:
                extra_kwargs["PyRosettaCluster_toml"] = self.toml
            if self.toml_format:
                extra_kwargs["PyRosettaCluster_toml_format"] = self.toml_format
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

    def _cache_toml(self) -> None:
        """Cache the Pixi/uv TOML file string and TOML file format."""

        if self.environment_manager == "pixi":
            # https://pixi.sh/dev/reference/environment_variables/#environment-variables-set-by-pixi
            toml_file = os.environ.get("PIXI_PROJECT_MANIFEST", "")
            if toml_file:
                if os.path.isfile(toml_file):
                    with open(toml_file, "r") as f:
                        self.toml = sanitize_urls(f.read())
                    self.toml_format = os.path.basename(toml_file)
                else:
                    logging.warning(
                        (
                            "`PyRosettaCluster` detected the set 'PIXI_PROJECT_MANIFEST' "
                            "environment variable, but the Pixi manifest file does not exist! "
                            "It is recommended to commit the Pixi manifest file to the "
                            "Git repository to reproduce the Pixi project later."
                        )
                    )
                    self.toml = ""
                    self.toml_format = ""
            else:
                # https://pixi.sh/dev/python/tutorial/#pixitoml-and-pyprojecttoml
                for filename in ("pixi.toml", "pyproject.toml"):
                    toml_file = os.path.join(os.getcwd(), filename)
                    if os.path.isfile(toml_file):
                        with open(toml_file, "r") as f:
                            self.toml = sanitize_urls(f.read())
                        self.toml_format = os.path.basename(toml_file)
                        break
                else:
                    logging.warning(
                        (
                            "`PyRosettaCluster` could not detect the Pixi manifest file! "
                            "It is recommended to commit the Pixi manifest file to the "
                            "Git repository to reproduce the Pixi project later."
                        )
                    )
                    self.toml = ""
                    self.toml_format = ""
        elif self.environment_manager == "uv":
            # https://docs.astral.sh/uv/reference/environment/#uv_project
            project_dir = os.environ.get("UV_PROJECT", None)
            if project_dir:
                toml_file = os.path.join(project_dir, "pyproject.toml")
                if os.path.isfile(toml_file):
                    with open(toml_file, "r") as f:
                        self.toml = sanitize_urls(f.read())
                    self.toml_format = os.path.basename(toml_file)
                else:
                    logging.warning(
                        (
                            "`PyRosettaCluster` detected the set 'UV_PROJECT' "
                            "environment variable, but the uv `pyproject.toml` file does not exist! "
                            "It is recommended to commit the uv `pyproject.toml` file to the "
                            "Git repository to reproduce the uv project later."
                        )
                    )
                    self.toml = ""
                    self.toml_format = ""
            else:
                toml_file = os.path.join(os.getcwd(), "pyproject.toml")
                if os.path.isfile(toml_file):
                    with open(toml_file, "r") as f:
                        self.toml = sanitize_urls(f.read())
                    self.toml_format = os.path.basename(toml_file)
                else:
                    logging.warning(
                        (
                            "`PyRosettaCluster` could not detect the uv `pyproject.toml` file! "
                            "The 'UV_PROJECT' environment variable is not set, and a `pyproject.toml` "
                            "file does not exist in the current working directory. For environment "
                            "reproducibility, please set the 'UV_PROJECT' environment variable "
                            "to the uv project root directory, or run the simulation from the uv project "
                            "root directory. If continuing with this simulation, it is recommended to commit the "
                            "uv `pyproject.toml` file to the Git repository to reproduce the uv project later."
                        )
                    )
                    self.toml = ""
                    self.toml_format = ""
        else:
            self.toml = ""
            self.toml_format = ""

    def _write_environment_file(self, filename: str) -> None:
        """
        Write the Conda/Mamba YML or uv/Pixi lock file string to the input filename. If Pixi/uv is used as the
        environment manager, also write the TOML file string to a separate filename.
        """

        if (
            (not self.simulation_records_in_scorefile)
            and (not self.dry_run)
            and self.environment
        ):
            with open(filename, "w") as f:
                f.write(self.environment)

            if self.environment_manager in ("pixi", "uv") and self.toml and self.toml_format:
                toml_file = "_".join(filename.split("_")[:-1] + [self.toml_format])
                with open(toml_file, "w") as f:
                    f.write(self.toml)

    def _write_init_file(self) -> None:
        """
        Maybe dump a PyRosetta initialization file.

        *Warning*: This method uses the `pickle` module to deserialize pickled `Pose` objects. Using the
        `pickle` module is not secure, so please only run with input files you trust. Learn more about the
        `pickle` module and its security `here <https://docs.python.org/3/library/pickle.html>`_.
        """

        if self.output_init_file != "":
            if self.compressed and self.output_init_file.endswith(".init"):
                self.output_init_file += ".bz2"
            if not self.output_init_file.endswith((".init", ".init.bz2")):
                raise ValueError(
                    "The output PyRosetta initialization file must end in '.init' or '.init.bz2'. "
                    + "The current `output_init_file` instance attribute of `PyRosettaCluster` is set to: "
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
        """
        Dump compressed PyRosetta initialization input files and `Pose` or `PackedPose` objects to the input
        filename.

        *Warning*: This method uses the `pickle` module to deserialize pickled `Pose` objects. Using the
        `pickle` module is not secure, so please only run with input files you trust. Learn more about the
        `pickle` module and its security `here <https://docs.python.org/3/library/pickle.html>`_.
        """

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
                "{0}: {1}. `{2}` did not run successfully, so PyRosetta initialization input files may not "
                "be saved! It is recommended to run `pyrosetta.dump_init_file` to reproduce PyRosetta "
                "initialization options on the head node later."
            )
            if self.compressed:
                try:
                    writer = PyRosettaInitDictWriter(**dump_init_file_kwargs)
                    init_file_json = writer.get_json()
                    if verbose:
                        writer.print_cached_files(filename, dump_init_file_kwargs["dry_run"])
                    with open(filename, "wb") as f:
                        f.write(bz2.compress(str.encode(init_file_json)))
                    logging.debug(_logging_debug_msg.format("PyRosettaInitDictWriter.get_json"))
                except Exception as ex:
                    logging.error(_logging_error_msg.format(type(ex).__name__, ex, "PyRosettaInitDictWriter.get_json"))
            else:
                try:
                    pyrosetta.dump_init_file(filename, **dump_init_file_kwargs)
                    logging.debug(_logging_debug_msg.format("pyrosetta.dump_init_file"))
                except Exception as ex:
                    logging.error(_logging_error_msg.format(type(ex).__name__, ex, "pyrosetta.dump_init_file"))
        out.flush()
        err.flush()


def verify_init_file(
    init_file: str,
    input_packed_pose: Optional[PackedPose],
    output_packed_pose: Optional[PackedPose],
    metadata: Dict[str, Any],
) -> None:
    """
    Verify that a PyRosetta initialization file was written by `PyRosettaCluster`.

    *Warning*: This function uses the `pickle` module to deserialize pickled `Pose` objects. Using the `pickle`
    module is not secure, so please only run with input files you trust. Learn more about the `pickle` module
    and its security `here <https://docs.python.org/3/library/pickle.html>`_.
    """

    @toolz.functoolz.curry
    def _verify_signer(_signer: InitFileSigner, _sha256: str, _signature: str) -> None:
        """Verify that current PyRosetta and `PyRosettaCluster` versions match that dumped in the '.init' file"""

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
    """
    Sign PyRosetta initialization file "metadata" and "poses" keys.

    *Warning*: This function uses the `pickle` module to deserialize pickled `Pose` objects. Using the `pickle`
    module is not secure, so please only run with input files you trust. Learn more about the `pickle` module
    and its security `here <https://docs.python.org/3/library/pickle.html>`_.
    """

    metadata = {}
    metadata["comment"] = "Generated by PyRosettaCluster"
    metadata["version"] = pyrosetta.distributed.cluster.__version__
    poses = []
    if isinstance(output_packed_pose, PackedPose): # Set output decoy to front of poses list if present
        _comments = dict(get_all_comments(output_packed_pose.pose))
        _key = IO.REMARK_FORMAT.rstrip() # Remove extra space since `add_comment` adds a space
        if _key not in _comments.keys():
            raise NotImplementedError(
                "Signing '.init' file metadata and poses when `PyRosettaCluster` simulation data is not "
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
) -> Tuple[Optional[PackedPose], Optional[PackedPose]]:
    """
    Return a `tuple` object of the input `PackedPose` object and the output `PackedPose` object from a ".init"
    or ".init.bz2" file, and optionally verify `PyRosettaCluster` metadata in the ".init" or ".init.bz2" file.

    *Warning*: This function uses the `pickle` module to deserialize pickled `Pose` objects. Using the `pickle`
    module is not secure, so please only run with input files you trust. Learn more about the `pickle` module
    and its security `here <https://docs.python.org/3/library/pickle.html>`_.
    """

    @toolz.functoolz.curry
    def _maybe_to_packed(
        key: str, poses: List[str], metadata: Dict[str, str]
    ) -> Optional[PackedPose]:
        assert isinstance(metadata, dict), (
            "The value of the 'metadata' key must be a `dict` object. "
            + f"Received: {type(metadata)}"
        )
        assert isinstance(poses, list), (
            "The value of the 'poses' key must be a `list` object. "
            + f"Received: {type(poses)}"
        )
        if key in metadata:
            idx = metadata[key]
            if isinstance(idx, int):
                return io.to_packed(io.to_pose(poses[idx]))
            else:
                raise TypeError(
                    f"The value of the '{key}' key in the value of the 'metadata' key must be an `int` object. "
                    + f"Received: {type(idx)}"
                )
        else:
            return None

    if not was_init_called():
        raise PyRosettaIsNotInitializedError(
            f"Please first initialize PyRosetta with the 'init_file' argument value: '{init_file}'"
        )
    if not isinstance(init_file, str):
        raise TypeError(f"The 'init_file' argument value must be a `str` object. Received: {type(init_file)}")
    if init_file.endswith(".init.bz2"):
        with open(init_file, "rb") as fbz2:
            init_dict = PyRosettaInitFileReader.from_json(bz2.decompress(fbz2.read()).decode())
    elif init_file.endswith(".init"):
        init_dict = PyRosettaInitFileReader.read_json(init_file)
    else:
        raise ValueError("The 'init_file' argument value must end with '.init' or '.init.bz2'.")

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
    Proxy for `pandas.read_pickle` for file-like objects using the `SecureSerializerBase` class in PyRosetta.
    Usage requires adding "pandas" as a secure package to unpickle in PyRosetta.

    *Warning*: This function uses the `pickle` module to deserialize pickled `pandas.DataFrame` objects. Using
    the `pickle` module is not secure, so please only run with input files you trust. Learn more about the
    `pickle` module and its security `here <https://docs.python.org/3/library/pickle.html>`_.

    Args:
        `filepath_or_buffer`: `str`
            See `pandas.read_pickle <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_pickle.html>`_.

        `compression`: `str | dict[str, Any] | None`
            See `pandas.read_pickle <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_pickle.html>`_.

            Default: `"infer"`

        `storage_options`: `dict[str, Any] | None`
            See `pandas.read_pickle <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_pickle.html>`_.

            Default: `None`

    Example:
        >>> pyrosetta.secure_unpickle.add_secure_package("pandas")
        >>> secure_read_pickle("/path/to/my/scorefile.gz")

    Note:
        If using `pandas` version `>=3.0.0`, PyArrow-backed datatypes may be enabled by default; in this case,
        please ensure that `pyrosetta.secure_unpickle.add_secure_package("pyarrow")` has also first been run.

        See https://pandas.pydata.org/pdeps/0010-required-pyarrow-dependency.html and
        https://pandas.pydata.org/pdeps/0014-string-dtype.html for more information.

    Returns:
        A deserialized `pandas.DataFrame` object.
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
    Scan the input string and sanitize any URLs that include credentials for source domains, returning the
    updated string.
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
                "`PyRosettaCluster` automatically removed embedded credentials from the "
                f"Conda channel '{host_domain}' while processing the environment file. "
                "These credentials are no longer required by this Conda channel. "
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


def _is_pandas_object_pyarrow_backed(obj: Union[pandas.DataFrame, pandas.Series]) -> bool:
    """
    Determine if a `pandas.DataFrame` or `pandas.Series` object uses Arrow-backed pandas dtypes.

    *Warning*: This function is experimental and subject to change in future versions.
    See https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.convert_dtypes.html
    for more information.

    Args:
        `obj`: `pandas.DataFrame | pandas.Series`
            An input `pandas.DataFrame` or `pandas.Series` object to test.

    Returns:
        A `bool` object.
    """

    def _is_arrow_dtype(dtype: Any) -> bool:
        return dtype.__class__.__name__ == "ArrowDtype"

    if isinstance(obj, pandas.DataFrame):
        if any(map(_is_arrow_dtype, obj.dtypes)):
            return True
        if any(map(_is_arrow_dtype, (obj.index.dtype, obj.columns.dtype))):
            return True
        return False
    elif isinstance(obj, pandas.Series):
        if any(map(_is_arrow_dtype, (obj.dtype, obj.index.dtype))):
            return True
        return False
    else:
        raise TypeError(f"Unsupported `pandas` object type: {type(obj)}")
