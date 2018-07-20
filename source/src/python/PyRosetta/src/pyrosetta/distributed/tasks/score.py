import threading

import traitlets
import pyrosetta.distributed
import pyrosetta.distributed.packed_pose as packed_pose
import pyrosetta.distributed.tasks.taskbase as taskbase
import pyrosetta.rosetta.core.scoring as scoring


class ScorePoseTask(taskbase.TaskBase):
    weights = traitlets.CUnicode(allow_none=True)
    patch = traitlets.CUnicode(allow_none=True)

    def __init__(self, weights=None, patch=None):
        super().__init__(weights=weights, patch=patch)

    def setup(self):
        self.protocol_lock = threading.Lock()

    @property
    @pyrosetta.distributed.requires_init
    @pyrosetta.distributed.with_lock
    def score_function(self):
        if not self.weights:
            return scoring.get_score_function()
        elif not self.patch:
            sff = scoring.ScoreFunctionFactory()
            return sff.create_score_function(self.weights)
        else:
            sff = scoring.ScoreFunctionFactory()
            return sff.create_score_function(
                self.weights, self.patch)

    @pyrosetta.distributed.requires_init
    def execute(self, pack_or_pose):
        work_pose = packed_pose.to_packed(pack_or_pose).pose
        with self.protocol_lock:
            self.score_function(work_pose)
        return packed_pose.PackedPose(work_pose)
