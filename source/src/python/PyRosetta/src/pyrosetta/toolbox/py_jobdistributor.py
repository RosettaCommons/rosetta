# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

import json
import os
import pyrosetta
import random
import warnings
from pyrosetta.bindings.pose import PoseScoreSerializer
from pyrosetta.rosetta.core.io.raw_data import ScoreMap
from sys import exit


_serializer = PoseScoreSerializer()

def output_scorefile(pose, pdb_name, current_name, scorefilepath, scorefxn, nstruct,
                     native_pose=None, additional_decoy_info=None, json_format=True):
    """
    Moved from PyJobDistributor (Jared Adolf-Bryfogle)
    Creates a scorefile if none exists, or appends the current one.
    Calculates and writes CA_rmsd if native pose is given,
    as well as any additional decoy info
    """


    if not json_format:
        exit(
            "\nThis scorefile output format is deprecated. Please use json output format instead or grab data from pose.scores")

    total_score = scorefxn(pose)  # Scores pose and calculates total score.
    entries = {}
    entries.update(
        {
         "pdb_name": str(pdb_name),
         "decoy": str(current_name),
         "filename": str(current_name),
         "nstruct": int(nstruct)
         }
    )
    scores =  dict(
        list(ScoreMap.get_arbitrary_string_data_from_pose(pose).items())
        + list(ScoreMap.get_arbitrary_score_data_from_pose(pose).items())
        + list(ScoreMap.get_energies_map_from_scored_pose(pose).items())
    )
    # Deserialize arbitrary python types from `ScoreMap.get_arbitrary_string_data_from_pose`
    scores = {scoretype: _serializer.maybe_decode(value) for scoretype, value in scores.items()}
    # Confirm that arbitrary python types can be JSON-encoded, otherwise remove them
    for key in list(scores.keys()):
        try:
            json.dumps(scores[key])
        except:
            warnings.warn(
                "Removing score key '{0}' with value of type '{1}' before ".format(key, type(scores[key]))
                + "saving the score file! Only JSON-serializable score values can be written to the "
                + "score file. Consider custom serializing the value to save this score or removing the "
                + "key from the `pose.scores` dictionary to remove this warning message."
            )
            scores.pop(key, None)

    entries.update(scores)
    if native_pose and "rmsd" not in entries:
        entries["rmsd"] = float(pyrosetta.rosetta.core.scoring.CA_rmsd(native_pose, pose))
    if additional_decoy_info:
        entries["additional_decoy_info"] = str(additional_decoy_info)

    with open(scorefilepath, "a") as f:
        json.dump(entries, f)
        f.write("\n")


class PyJobDistributor:

    def __init__(self, pdb_name, nstruct, scorefxn, compress=False):
        self.pdb_name = pdb_name
        self.nstruct = nstruct
        self.compress = compress

        self.current_id = None
        self.current_name = None          # Current decoy name
        self.scorefxn = scorefxn          # Used for final score calculation
        self.native_pose = None           # Used for rmsd calculation
        self.additional_decoy_info = None # Used for any additional decoy information you want stored
        self.json_format = True          # Used for JSON formatted scorefile

        self.sequence = list(range(nstruct))
        random.shuffle(self.sequence)
        self.start_decoy()                # Initializes the job distributor

    @property
    def job_complete(self):
        return len(self.sequence) == 0

    def start_decoy(self):
        while self.sequence:

            self.current_id = self.sequence[0]
            self.current_name = self.pdb_name + "_" + str(self.current_id) + (".pdb.gz" if self.compress else ".pdb")
            self.current_in_progress_name = self.current_name + ".in_progress"

            if (not os.path.isfile(self.current_name)) and (not os.path.isfile(self.current_in_progress_name)):
                with open(self.current_in_progress_name, "w") as f:
                    f.write("This decoy is in progress.")
                print("Working on decoy: {}".format(self.current_name))
                break

    def output_decoy(self, pose):
        if os.path.isfile(self.current_name): # decoy already exists, probably written to other process -> moving to next decoy if any

            if os.path.isfile(self.current_in_progress_name):
                os.remove(self.current_in_progress_name)

            if not self.job_complete:
                self.start_decoy()
                self.output_decoy(self, pose)

        else:
            if self.compress:
                s = pyrosetta.rosetta.std.ostringstream()
                pose.dump_pdb(s)

                z = pyrosetta.rosetta.utility.io.ozstream(self.current_name)
                z.write(s.str(), len(s.str()))
                del z

            else:
                pose.dump_pdb(self.current_name)

            score_tag = ".fasc" if pose.is_fullatom() else ".sc"

            scorefile = self.pdb_name + score_tag
            output_scorefile(pose, self.pdb_name, self.current_name, scorefile, self.scorefxn,
                             self.nstruct, self.native_pose, self.additional_decoy_info, self.json_format)

            self.sequence.remove(self.current_id)

            if os.path.isfile(self.current_in_progress_name):
                os.remove(self.current_in_progress_name)

            self.start_decoy()
