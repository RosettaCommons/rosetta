# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import types
import warnings

try:
    from collections.abc import MutableMapping
except ImportError:
    # For python < 3.3
    from collections import MutableMapping

from pyrosetta.rosetta.core.pose import (
    getPoseExtraFloatScores,
    getPoseExtraStringScores,
    setPoseExtraScore,
    hasPoseExtraScore,
    hasPoseExtraScore_str,
    clearPoseExtraScore,
    clearPoseExtraScores,
)

from pyrosetta.bindings.scores.base import PoseCacheAccessorBase


class ExtraScoresAccessorBase(PoseCacheAccessorBase, MutableMapping):
    """
    Base methods for accessor wrapper for pose arbitrary extra string scores
    and pose arbitrary extra float scores.
    """
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    def __getitem__(self, key):
        return self.maybe_decode(self.all[key])

    def __setitem__(self, key, value):
        self._validate_set(key)
        value = self.maybe_encode(value)

        if isinstance(value, str):
            if hasPoseExtraScore(self.pose, key):
                # User is setting an arbitrary string score, but key exists in arbitrary extra float scores
                raise KeyError(
                    "Cannot set an identical key that is being used "
                    + "in the arbitrary extra float scores: '{0}'".format(key)
                )
        elif isinstance(value, float):
            if hasPoseExtraScore_str(self.pose, key):
                # User is setting an arbitrary float score, but key exists in arbitrary extra string scores
                raise KeyError(
                    "Cannot set an identical key that is being used "
                    + "in the arbitrary extra string scores: '{0}'".format(key)
                )
        else:
            raise NotImplementedError(type(value))

        setPoseExtraScore(self.pose, key, value)

    def __delitem__(self, key):
        self._validate_del(key)
        if not hasPoseExtraScore(self.pose, key) and not hasPoseExtraScore_str(self.pose, key):
            raise KeyError(
                "Key is not in arbitrary extra string scores or arbitrary extra float scores: {0}".format(key)
            )
        clearPoseExtraScore(self.pose, key)

    def clear(self):
        for key in self.all.keys():
            clearPoseExtraScore(self.pose, key)


class ExtraFloatScoresDataAccessor(ExtraScoresAccessorBase):
    """Accessor wrapper for pose arbitrary extra float scores."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        return types.MappingProxyType(
            dict(getPoseExtraFloatScores(self.pose))
        )

    def __setitem__(self, key, value):
        if not isinstance(value, float):
            warnings.warn(
                "Key '{0}' with value of type '{1}' is being saved as arbitrary extra string data.".format(
                    key, type(value),
                ),
                UserWarning,
                stacklevel=2,
            )
        return super().__setitem__(key, value)


class ExtraStringScoresDataAccessor(ExtraScoresAccessorBase):
    """Accessor wrapper for pose arbitrary extra string scores."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        return types.MappingProxyType(
            dict(getPoseExtraStringScores(self.pose))
        )

    def __setitem__(self, key, value):
        if isinstance(value, float):
            warnings.warn(
                "Key '{0}' with value of type '{1}' is being saved as arbitrary extra real data.".format(
                    key, type(value),
                ),
                UserWarning,
                stacklevel=2,
            )
        return super().__setitem__(key, value)


class ExtraScoresAccessor(ExtraScoresAccessorBase):
    """Accessor wrapper for pose arbitrary extra string scores and pose arbitrary extra float scores."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        """
        Get all arbitrary extra float and extra string scores (with clobber warnings).

        This method aims to mimic data override precedences used in the legacy `pose.scores` dictionary:
            1. `pose.energies().active_total_energies()`
            2. `ScoreMap.get_arbitrary_score_data_from_pose(pose)`
            3. `ScoreMap.get_arbitrary_string_data_from_pose(pose)`
        """
        extra_string_scores = self.string
        extra_float_scores = self.real

        for k in extra_float_scores.keys():
            if k in extra_string_scores.keys():
                self._clobber_warning(
                    "Arbitrary extra float score key is clobbering arbitrary extra string score key: '{0}'".format(k)
                )

        return types.MappingProxyType(
            {
                **extra_string_scores,
                **extra_float_scores,
            }
        )

    @property
    def string(self):
        return ExtraStringScoresDataAccessor(self.pose)

    @property
    def real(self):
        return ExtraFloatScoresDataAccessor(self.pose)

    def clear(self):
        """Clear pose extra scores."""
        clearPoseExtraScores(self.pose)
