# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


import warnings

from pyrosetta.rosetta.core.simple_metrics import clear_sm_data, get_sm_data
from pyrosetta.rosetta.core.simple_metrics.metrics import (
    CustomRealValueMetric,
    CustomStringValueMetric,
)
from pyrosetta.rosetta.core.scoring import ScoreType

from pyrosetta.bindings.scores.serialization import PoseScoreSerializer


class ClobberWarning(UserWarning):
    pass


class PoseCacheAccessorBase(PoseScoreSerializer):
    __slots__ = ("pose",)

    def __init__(self, pose):
        self.pose = pose
        self.custom_real_value_metric = CustomRealValueMetric()
        self.custom_string_value_metric = CustomStringValueMetric()

    @property
    def _reserved_custom_metric_keys(self):
        return {"custom_real_valued_metric", "custom_string_valued_metric"}

    @property
    def _reserved(self):
        _reserved = {k for k in ScoreType.__dict__.keys() if not k.startswith("_")}
        _reserved = _reserved.union(self._reserved_custom_metric_keys)
        return _reserved

    @property
    def _sm_data_accessor_attrs(self):
        """Supported `SimpleMetricDataAccessor` attributes to reset after `clear_sm_data`."""
        return (
            "string",
            "real",
        )

    @property
    def _unsupported_sm_data_accessor_attrs(self):
        """Unsupported `SimpleMetricDataAccessor` attributes that cannot be reset after `clear_sm_data`."""
        return (
            "composite_string",
            "composite_real",
            "per_residue_string",
            "per_residue_real",
            "per_residue_probabilities",
        )

    def _get_sm_data_dict(self, attributes):
        return {
            _attr: dict(getattr(self.pose.cache.metrics, _attr))
            for _attr in attributes
        }

    def _maybe_delete_keys_from_sm_data(self, keys=None, attributes=None):
        """
        Cache, clear, and restore all SimpleMetric data except the provided keys
        in the provided `SimpleMetricDataAccessor` attributes. This is necessary to
        delete one or more keys, since `clear_sm_data` clears all SimpleMetric data.
        """
        keys = tuple() if keys is None else tuple(keys)
        attributes = tuple() if attributes is None else tuple(attributes)
        if any(_attr not in self._sm_data_accessor_attrs for _attr in attributes):
            raise NotImplementedError(attributes)
        sm_data = get_sm_data(self.pose).get_all_sm_data()
        if sm_data.has_data():
            # Cache SimpleMetric data
            _sm_data_dict = self._get_sm_data_dict(self._sm_data_accessor_attrs)
            _unsupported_sm_data_dict = self._get_sm_data_dict(self._unsupported_sm_data_accessor_attrs)
            # Raise if deleting an item we cannot restore
            for _d in _unsupported_sm_data_dict.values():
                for _k in _d.keys():
                    for key in keys:
                        if key == _k:
                            raise KeyError(
                                "Cannot delete a SimpleMetric data key '{0}' from: {1}. ".format(
                                    key, self._unsupported_sm_data_accessor_attrs
                                ) + "Consider using `pose.cache.metrics.clear()`."
                            )
            if any(map(len, _unsupported_sm_data_dict.values())):
                raise KeyError(
                    "Cannot delete SimpleMetric data keys {0} because Pose contains SimpleMetric data from: {1}".format(
                        keys, self._unsupported_sm_data_accessor_attrs
                    ) + "Consider using `pose.cache.metrics.clear()`."
                )
            # Clear SimpleMetric data
            clear_sm_data(self.pose)
            # Reapply all SimpleMetric data except the deleted item
            for _attr, _d in _sm_data_dict.items():
                for _k, _v in _d.items():
                    if not (_k in keys and _attr in attributes):
                        if _attr == "string":
                            m = self.custom_string_value_metric
                        elif _attr == "real":
                            m = self.custom_real_value_metric
                        else:
                            raise NotImplementedError(_attr)
                        m.set_value(self.maybe_encode(_v))
                        m.apply(out_label=_k, pose=self.pose, override_existing_data=True)

    def _maybe_delete_reserved_keys_from_sm_data(self):
        """
        `CustomRealValueMetric` and `CustomStringValueMetric` can save the extra default keys
        "custom_real_valued_metric" or "custom_string_valued_metric" because `pose.cache.metrics`
        does not set a 'custom_type' parameter when setting SimpleMetric data. This method
        aims to delete these extra default keys if they were created and we can delete them, 
        otherwise fallback to warning the user that they were created.
        """
        _has_custom_string_valued_metric = "custom_string_valued_metric" in self.pose.cache.metrics.string.keys()
        _has_custom_real_valued_metric = "custom_real_valued_metric" in self.pose.cache.metrics.real.keys()
        _static_keys = (
            self.pose.cache.metrics.composite_string.keys(),
            self.pose.cache.metrics.composite_real.keys(),
            self.pose.cache.metrics.per_residue_string.keys(),
            self.pose.cache.metrics.per_residue_real.keys(),
            self.pose.cache.metrics.per_residue_probabilities.keys(),
        )
        if not any(map(len, _static_keys)):
            if _has_custom_string_valued_metric and _has_custom_real_valued_metric:
                self._maybe_delete_keys_from_sm_data(
                    keys=("custom_string_valued_metric", "custom_real_valued_metric",),
                    attributes=("string", "real"),
                )
            elif _has_custom_string_valued_metric:
                self._maybe_delete_keys_from_sm_data(
                    keys=("custom_string_valued_metric",),
                    attributes=("string",),
                )
            elif _has_custom_real_valued_metric:
                self._maybe_delete_keys_from_sm_data(
                    keys=("custom_real_valued_metric",),
                    attributes=("real",)
                )
        else:
            if _has_custom_string_valued_metric:
                warnings.warn(
                    "`CustomStringValueMetric` implemented the key 'custom_string_valued_metric' in `pose.cache`.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            if _has_custom_real_valued_metric:
                warnings.warn(
                    "`CustomRealValueMetric` implemented the key 'custom_real_valued_metric' in `pose.cache`.",
                    RuntimeWarning,
                    stacklevel=2,
                )

    def _clobber_warning(self, msg):
        warnings.warn(
            msg,
            ClobberWarning,
            stacklevel=2,
        )

    def _validate_set(self, key):
        if key in self._reserved:
            raise KeyError("Cannot set a key with a reserved name: {0}".format(key))

    def _validate_del(self, key):
        if key in self._reserved:
            raise KeyError("Cannot delete a key with a reserved name: {0}".format(key))

    def __len__(self):
        return len(self.all)

    def __iter__(self):
        return iter(self.all)

    def __str__(self):
        return str(dict(self))

    def _repr_pretty_(self, p, cycle):
        """IPython-display representation."""
        p.pretty(dict(self))
