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
import tempfile
import uuid

from contextlib import redirect_stdout, redirect_stderr
from datetime import datetime
from pyrosetta.rosetta.core.pose import Pose, add_comment
from pyrosetta.distributed.packed_pose.core import PackedPose
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

from pyrosetta.distributed.cluster.exceptions import OutputError
from pyrosetta.distributed.cluster.init_files import setup_init_file_metadata_and_poses
from pyrosetta.distributed.cluster.serialization import update_scores
from pyrosetta.distributed.cluster.logging_support import RedirectToLogger


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
                json.dumps(scores_dict[key])
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
            pdbfile_data = json.dumps(simulation_data_json, sort_keys=False, indent=None, separators=(", ", ": "))
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
                output_init_file = os.path.join(output_dir, decoy_name + ".init")
                if self.compressed:
                    with tempfile.TemporaryDirectory() as tmp_dir:
                        _tmp_init_file = os.path.join(tmp_dir, "tmp.init")
                        self._dump_init_file(
                            _tmp_init_file,
                            input_packed_pose=self.input_packed_pose,
                            output_packed_pose=_packed_pose,
                            verbose=False,
                        )
                        with open(_tmp_init_file, "r") as f:
                            init_file_str = f.read()
                    output_init_file += ".bz2"
                    with open(output_init_file, "wb") as f:
                        f.write(bz2.compress(str.encode(init_file_str)))
                else:
                    self._dump_init_file(
                        output_init_file,
                        input_packed_pose=self.input_packed_pose,
                        output_packed_pose=_packed_pose,
                        verbose=False,
                    )

            # Output JSON-encoded scorefile
            if ".json" in self.output_scorefile_types:
                if self.simulation_records_in_scorefile:
                    scorefile_data = pdbfile_data
                else:
                    scorefile_data = json.dumps(
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
                        df_chunk = pandas.read_pickle(_scorefile_path, compression="infer")
                        df = pandas.concat([df_chunk, df])
                    df.to_pickle(_scorefile_path, compression="infer")

    def _write_environment_file(self, filename: str) -> None:
        """Write the YML string to the input filename."""

        if (
            (not self.simulation_records_in_scorefile)
            and (not self.dry_run)
            and self.environment
        ):
            with open(filename, "w") as f:
                f.write(self.environment)

    def _write_init_file(self, filename: str) -> None:
        """Maybe write PyRosetta initialization file to the input filename."""

        if filename != "":
            self._dump_init_file(
                filename,
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
            metadata, poses = setup_init_file_metadata_and_poses(
                input_packed_pose=input_packed_pose,
                output_packed_pose=output_packed_pose,
            )
            try:
                pyrosetta.dump_init_file(
                    filename,
                    poses=poses,
                    author=self.author,
                    email=self.email,
                    license=self.license,
                    metadata=metadata,
                    overwrite=True,
                    dry_run=self.dry_run,
                    verbose=verbose,
                )
                logging.debug("Successfully ran `pyrosetta.dump_init_file` on the host node.")
            except Exception as ex:
                logging.error(
                    f"{type(ex).__name__}: {ex}. `pyrosetta.dump_init_file` did not run successfully, "
                    + "so PyRosetta initialization input data may not be saved! It is recommended to run "
                    + "`pyrosetta.dump_init_file` to reproduce PyRosetta initialization options on the "
                    + "host node later."
                )
        out.flush()
        err.flush()
