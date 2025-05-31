# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

try:
    import toolz
except ImportError:
    print(
        "Importing 'pyrosetta.distributed.cluster.io' requires the "
        + "third-party package 'toolz' as a dependency!\n"
        + "Please install this package into your python environment. "
        + "For installation instructions, visit:\n"
        + "https://pypi.org/project/toolz/\n"
    )
    raise

import bz2
import collections
import json
import logging
import os
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import uuid

from datetime import datetime
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.distributed.cluster.exceptions import OutputError
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
        instance_kwargs = self.serializer.deepcopy_kwargs(instance_state)
        for i in self.__attrs_attrs__:
            if not i.init:
                instance_kwargs.pop(i.name)
            if i.name == "input_packed_pose":
                instance_kwargs.pop(i.name, None)
        instance_kwargs["tasks"] = kwargs.pop("task")
        for option in ["extra_options", "options"]:
            if option in instance_kwargs["tasks"]:
                instance_kwargs["tasks"][option] = pyrosetta.distributed._normflags(
                    instance_kwargs["tasks"][option]
                )
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
                    + "key from the `pose.scores` dictionary to remove this warning message."
                )
                scores_dict.pop(key, None)

        return scores_dict

    @staticmethod
    def _format_result(result: Union[Pose, PackedPose]) -> Tuple[str, Dict[Any, Any]]:
        """
        Given a `Pose` or `PackedPose` object, return a tuple containing
        the pdb string and a scores dictionary.
        """

        _pdbstring = io.to_pdbstring(result)
        _scores_dict = io.to_dict(result)
        _scores_dict.pop("pickled_pose", None)
        _scores_dict = IO._filter_scores_dict(_scores_dict)

        return (_pdbstring, _scores_dict)

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
                out = [IO._format_result(results)]
            else:
                out = []
        elif isinstance(results, collections.abc.Iterable):
            out = []
            for result in results:
                if isinstance(results, bytes):
                    result = self.serializer.decompress_packed_pose(result)
                if isinstance(result, (Pose, PackedPose)):
                    if not io.to_pose(result).empty():
                        out.append(IO._format_result(result))
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

    def _save_results(self, results: Any, kwargs: Dict[Any, Any]) -> None:
        """Write results and kwargs to disk."""

        if self.dry_run:
            logging.info(
                "PyRosettaCluster `dry_run` attribute is enabled. Skipping saving!"
            )
            return

        # Parse and save results
        for pdbstring, scores in self._parse_results(results):
            kwargs = self._process_kwargs(kwargs)
            output_dir = self._get_output_dir(decoy_dir=self.decoy_path)
            decoy_name = "_".join([self.simulation_name, uuid.uuid4().hex])
            output_file = os.path.join(output_dir, decoy_name + ".pdb")
            if self.compressed:
                output_file += ".bz2"
            extra_kwargs = {
                "PyRosettaCluster_decoy_name": decoy_name,
                "PyRosettaCluster_output_file": output_file,
            }
            if os.path.exists(self.environment_file):
                extra_kwargs[
                    "PyRosettaCluster_environment_file"
                ] = self.environment_file
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
            pdbfile_data = json.dumps(
                {
                    "instance": collections.OrderedDict(sorted(instance.items())),
                    "metadata": collections.OrderedDict(sorted(metadata.items())),
                    "scores": collections.OrderedDict(sorted(scores.items())),
                }
            )
            # Write full .pdb record
            pdbstring_data = pdbstring + os.linesep + self.REMARK_FORMAT + pdbfile_data
            if self.compressed:
                with open(output_file, "wb") as f:
                    f.write(bz2.compress(str.encode(pdbstring_data)))
            else:
                with open(output_file, "w") as f:
                    f.write(pdbstring_data)
            if self.simulation_records_in_scorefile:
                scorefile_data = pdbfile_data
            else:
                scorefile_data = json.dumps(
                    {
                        metadata["output_file"]: collections.OrderedDict(
                            sorted(scores.items())
                        ),
                    }
                )
            # Append data to scorefile
            with open(self.scorefile_path, "a") as f:
                f.write(scorefile_data + os.linesep)

    def _write_environment_file(self, filename: str) -> None:
        """Write the YML string to the input filename."""

        if (
            (not self.simulation_records_in_scorefile)
            and (not self.dry_run)
            and self.environment
        ):
            with open(filename, "w") as f:
                f.write(self.environment)
