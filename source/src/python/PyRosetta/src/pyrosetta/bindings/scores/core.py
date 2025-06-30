# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import types

try:
    from collections.abc import MutableMapping
except ImportError:
    # For python < 3.3
    from collections import MutableMapping

from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.pose import (
    getPoseExtraFloatScores,
    getPoseExtraStringScores,
    setPoseExtraScore,
    hasPoseExtraScore,
    hasPoseExtraScore_str,
    clearPoseExtraScore,
    clearPoseExtraScores,
)
from pyrosetta.rosetta.core.simple_metrics import clear_sm_data, get_sm_data

from pyrosetta.bindings.scores.base import PoseCacheAccessorBase
from pyrosetta.bindings.scores.extra_scores import ExtraScoresAccessor
from pyrosetta.bindings.scores.simple_metrics import SimpleMetricDataAccessor
from pyrosetta.bindings.scores.energies import EnergiesAccessor


class PoseCacheAccessor(PoseCacheAccessorBase, MutableMapping):
    """Accessor wrapper for pose energies, extra scores, and SimpleMetrics."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def extra(self):
        return ExtraScoresAccessor(self.pose)

    @property
    def metrics(self):
        return SimpleMetricDataAccessor(self.pose)

    @property
    def energies(self):
        return EnergiesAccessor(self.pose)

    @property
    def all_scores(self):
        """Return a nested dictionary of all score data."""
        return types.MappingProxyType(
            {
                "extra": types.MappingProxyType(
                    {
                        "string": types.MappingProxyType(dict(self.extra.string)),
                        "real": types.MappingProxyType(dict(self.extra.real)),
                    }
                ),
                "metrics": types.MappingProxyType(
                        {
                        "string": types.MappingProxyType(dict(self.metrics.string)),
                        "real": types.MappingProxyType(dict(self.metrics.real)),
                        "composite_string": types.MappingProxyType(dict(self.metrics.composite_string)),
                        "composite_real": types.MappingProxyType(dict(self.metrics.composite_real)),
                        "per_residue_string": types.MappingProxyType(dict(self.metrics.per_residue_string)),
                        "per_residue_real": types.MappingProxyType(dict(self.metrics.per_residue_real)),
                        "per_residue_probabilities": types.MappingProxyType(dict(self.metrics.per_residue_probabilities)),
                    }
                ),
                "energies": types.MappingProxyType(dict(self.energies)),
            }
        )

    @property
    def all_keys(self):
        """Return a tuple of all score data keys."""
        ref_keys = []
        for k1, v1 in self.all_scores.items():
            for k2, v2 in v1.items():
                if k1 == "energies":
                    ref_keys.append(k2)
                else:
                    for k3 in v2.keys():
                        ref_keys.append(k3)

        return tuple(ref_keys)

    def assert_unique_keys(self):
        """Assert that `Pose.cache.all_keys` returns all unique keys."""
        ref_keys = self.all_keys
        if len(ref_keys) == len(set(ref_keys)):
            print("Pose cached data keys are unique: {0}".format(ref_keys))
        else:
            found = set()
            duplicate_keys = {k for k in ref_keys if k in found or found.add(k)}
            raise AssertionError(
                "Pose cached data keys are not unique: {0}".format(duplicate_keys)
            )

    @property
    def all(self):
        """
        This method aims to mimic data override precedences used in the legacy `pose.scores` dictionary:
            1. `pose.energies().active_total_energies()`
            2. `ScoreMap.get_arbitrary_score_data_from_pose(pose)`
            3. `ScoreMap.get_arbitrary_string_data_from_pose(pose)`
        Data override precedences as defined in `ScoreMap::add_arbitrary_score_data_from_pose`:
            1. SimpleMetric data
            2. Arbitrary extra float scores
        Data override precedences as defined in `ScoreMap::add_arbitrary_string_data_from_pose`:
            1. SimpleMetric data
            2. Arbitrary extra string scores
        """
        extra = self.extra
        metrics = self.metrics
        energies = self.energies

        for k in extra.keys():
            if k in metrics.keys():
                self._clobber_warning(
                    "SimpleMetrics data key is clobbering arbitrary extra scores data key: '{0}'".format(k)
                )
            if k in energies.keys():
                self._clobber_warning(
                    "Active total energies key is clobbering arbitrary extra scores data key: '{0}'".format(k)
                )
        for k in metrics.keys():
            if k in energies.keys():
                self._clobber_warning(
                    "Active total energies key is clobbering SimpleMetrics data key: '{0}'".format(k)
                )

        return types.MappingProxyType(
            {
                **extra,
                **metrics,
                **energies,
            }
        )

    def __getitem__(self, key):
        "Get a value from a key from the `Pose.cache` dictionary."
        return self.maybe_decode(self.all[key])

    def __setitem__(self, key, value):
        """
        Save a key and value to SimpleMetric data by default.

        To save a key and value to arbitrary extra scores, use:
            `Pose.cache.extra[key] = value`
        """
        self.metrics.__setitem__(key, value)

    def __delitem__(self, key):
        """
        Delete an item from the pose scores cache. This method first looks for the key in SimpleMetric
        data, then if not present it looks for the key in arbitrary extra float data and arbitrary extra
        string data. Use `pose.cache.extra.pop(key)` to explicitly delete an item from the arbitrary extra
        float data or arbitrary extra string data. Use `pose.cache.metrics.pop(key)` to explicitly delete
        an item from SimpleMetric data.
        """
        self._validate_del(key)
        sm_data = get_sm_data(self.pose).get_all_sm_data()
        if sm_data.has_data() and key in self.metrics.keys():
            self._maybe_delete_keys_from_sm_data(keys=(key,), attributes=("string", "real"))
        elif hasPoseExtraScore(self.pose, key) or hasPoseExtraScore_str(self.pose, key):
            clearPoseExtraScore(self.pose, key)
        else:
            raise KeyError(key)

    def clear(self):
        """Clear pose energies, extra scores, and SimpleMetric data."""
        self.pose.energies().clear()
        clearPoseExtraScores(self.pose)
        clear_sm_data(self.pose)
