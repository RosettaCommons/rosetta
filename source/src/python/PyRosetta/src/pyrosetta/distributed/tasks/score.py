import threading

import traitlets
import pyrosetta.distributed
import pyrosetta.distributed.packed_pose as packed_pose
import pyrosetta.distributed.tasks.taskbase as taskbase
import pyrosetta.rosetta.core.scoring as scoring

from pyrosetta.rosetta.protocols.loops import get_fa_scorefxn


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
            sff = scoring.ScoreFunctionFactory.get_instance()
            return sff.create_score_function(self.weights)
        else:
            sff = scoring.ScoreFunctionFactory.get_instance()
            return sff.create_score_function(
                self.weights, self.patch)

    @pyrosetta.distributed.requires_init
    def execute(self, pack_or_pose):
        work_pose = packed_pose.to_packed(pack_or_pose).pose
        with self.protocol_lock:
            self.score_function(work_pose)
        return packed_pose.PackedPose(work_pose)


def create_score_function(*args, **kwargs):
    """
    Returns a `ScorePoseTask` instance.
    Input `*args` and `**kwargs` are passed to `ScorePoseTask`.

    @klimaj
    """
    return ScorePoseTask(*args, **kwargs)


def get_fa_scorefxn():
    """
    Returns a `ScorePoseTask` instance with `weights`
    from `pyrosetta.rosetta.protocols.loops.get_fa_scorefxn()`.

    @klimaj
    """
    weights = get_fa_scorefxn().get_name()
    return ScorePoseTask(weights=weights, patch=None)


def get_score_function():
    """
    Returns a `ScorePoseTask` instance with `weights` and `patch`
    from `pyrosetta.rosetta.core.scoring.get_score_function()`.

    @klimaj
    """
    return ScorePoseTask(weights=None, patch=None)
