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

from pyrosetta.rosetta.core.simple_metrics import clear_sm_data, get_sm_data

from pyrosetta.bindings.scores.base import PoseCacheAccessorBase


class SimpleMetricDataAccessorBase(PoseCacheAccessorBase, MutableMapping):
    """Base methods for accessor wrapper for pose SimpleMetric data."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    def __getitem__(self, key):
        return self.maybe_decode(self.all[key])

    def __setitem__(self, key, value):
        self._validate_set(key)
        if isinstance(value, float):
            m = self.custom_real_value_metric
        else:
            m = self.custom_string_value_metric
        m.set_value(self.maybe_encode(value))
        m.apply(out_label=key, pose=self.pose, override_existing_data=True)
        self._maybe_delete_reserved_keys_from_sm_data()

    def __delitem__(self, key):
        self._validate_del(key)
        sm_data = get_sm_data(self.pose).get_all_sm_data()
        if sm_data.has_data():
            self._maybe_delete_keys_from_sm_data(keys=(key,), attributes=("string", "real"))
        else:
            raise KeyError(key)

    def _format_metric(self, raw_data, as_dict=False):
        out_data = {} 
        for name, data in raw_data.items():
            for first, second in data.items():
                key = "_".join([name, first])
                out_data[key] = dict(second) if as_dict else second
        return out_data

    def format_composite_string(self, raw_data):
        """Mimics formatting in `ScoreMap.add_arbitrary_string_data_from_pose`."""
        return self._format_metric(raw_data, as_dict=False)

    def format_composite_real(self, raw_data):
        """Mimics formatting in `ScoreMap.add_arbitrary_score_data_from_pose`."""
        return self._format_metric(raw_data, as_dict=False)

    def format_per_residue_string(self, raw_data):
        """Mimics formatting in `ScoreMap.add_arbitrary_string_data_from_pose`."""
        return self._format_metric(raw_data, as_dict=False)

    def format_per_residue_real(self, raw_data):
        """Mimics formatting in `ScoreMap.add_arbitrary_score_data_from_pose`."""
        return self._format_metric(raw_data, as_dict=False)

    def format_per_residue_probabilities(self, raw_data):
        """
        Format per-residue probabilities metrics as nested dictionaries.
        
        The default output format from `ScoreMap.add_arbitrary_string_data_from_pose`
        sets values as `str` objects:
            {'my_metric_1': 'ALA:0.0125,ASN:0.000000,...', 'my_metric_2': ...}`

        Here, we reformat the values into separate `dict` objects where keys are `str`
        objects and values are `float` objects:
            {'my_metric_1': {'ALA': 0.0125, 'ASN': 0.0, ...}, 'my_metric_2': ...}
        """
        return self._format_metric(raw_data, as_dict=True)

    def clear(self):
        self._maybe_delete_keys_from_sm_data(
            keys=tuple(self.all.keys()),
            attributes=("string", "real"),
        )


class SimpleMetricStringDataAccessor(SimpleMetricDataAccessorBase):
    """Accessor wrapper for pose SimpleMetric string data."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        return types.MappingProxyType(
            dict(get_sm_data(self.pose).get_string_metric_data())
        )

    def __setitem__(self, key, value):
        if isinstance(value, float):
            warnings.warn(
                "Key '{0}' with value of type '{1}' is being saved as SimpleMetric real data.".format(
                    key, type(value),
                ),
                UserWarning,
                stacklevel=2,
            )
        return super().__setitem__(key, value)

    def __delitem__(self, key):
        self._validate_del(key)
        sm_data = get_sm_data(self.pose).get_all_sm_data()
        if sm_data.has_data():
            self._maybe_delete_keys_from_sm_data(keys=(key,), attributes=("string",))
        else:
            raise KeyError(key)

    def clear(self):
        self._maybe_delete_keys_from_sm_data(keys=tuple(self.all.keys()), attributes=("string",))


class SimpleMetricRealDataAccessor(SimpleMetricDataAccessorBase):
    """Accessor wrapper for pose SimpleMetric real data."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        return types.MappingProxyType(
            dict(get_sm_data(self.pose).get_real_metric_data())
        )

    def __setitem__(self, key, value):
        if not isinstance(value, float):
            warnings.warn(
                "Key '{0}' with value of type '{1}' is being saved as SimpleMetric string data.".format(
                    key, type(value),
                ),
                UserWarning,
                stacklevel=2,
            )
        return super().__setitem__(key, value)

    def __delitem__(self, key):
        self._validate_del(key)
        sm_data = get_sm_data(self.pose).get_all_sm_data()
        if sm_data.has_data():
            self._maybe_delete_keys_from_sm_data(keys=(key,), attributes=("real",))
        else:
            raise KeyError(key)

    def clear(self):
        self._maybe_delete_keys_from_sm_data(keys=tuple(self.all.keys()), attributes=("real",))


class SimpleMetricCompositeStringDataAccessor(SimpleMetricDataAccessorBase):
    """Accessor wrapper for pose SimpleMetric composite string data."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        return types.MappingProxyType(
            self.format_composite_string(
                get_sm_data(self.pose).get_composite_string_metric_data()
            )
        )

    def __setitem__(self, key, value):
        raise NotImplementedError("Cannot set SimpleMetric composite string data.")


class SimpleMetricCompositeRealDataAccessor(SimpleMetricDataAccessorBase):
    """Accessor wrapper for pose SimpleMetric composite real data."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        return types.MappingProxyType(
            self.format_composite_real(
                get_sm_data(self.pose).get_composite_real_metric_data()
            )
        )

    def __setitem__(self, key, value):
        raise NotImplementedError("Cannot set SimpleMetric composite real data.")


class SimpleMetricPerResidueStringDataAccessor(SimpleMetricDataAccessorBase):
    """Accessor wrapper for pose SimpleMetric per-residue string data."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        return types.MappingProxyType(
            self.format_per_residue_string(
                get_sm_data(self.pose).get_per_residue_string_metric_output()
            )
        )

    def __setitem__(self, key, value):
        raise NotImplementedError("Cannot set SimpleMetric per-residue string data.")


class SimpleMetricPerResidueRealDataAccessor(SimpleMetricDataAccessorBase):
    """Accessor wrapper for pose SimpleMetric per-residue real data."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        return types.MappingProxyType(
            self.format_per_residue_real(
                get_sm_data(self.pose).get_per_residue_real_metric_output()
            )
        )

    def __setitem__(self, key, value):
        raise NotImplementedError("Cannot set SimpleMetric per-residue real data.")


class SimpleMetricPerResidueProbabilitiesDataAccessor(SimpleMetricDataAccessorBase):
    """Accessor wrapper for pose SimpleMetric per-residue probabilities data."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        return types.MappingProxyType(
            self.format_per_residue_probabilities(
                get_sm_data(self.pose).get_per_residue_probabilities_metric_output()
            )
        )

    def __setitem__(self, key, value):
        raise NotImplementedError("Cannot set SimpleMetric per-residue probabilities data.")


class SimpleMetricDataAccessor(SimpleMetricDataAccessorBase):
    """Accessor wrapper for pose SimpleMetric data."""
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        """
        This method aims to mimic data override precedences used in the legacy `pose.scores` dictionary:
            1. `pose.energies().active_total_energies()`
            2. `ScoreMap.get_arbitrary_score_data_from_pose(pose)`
            3. `ScoreMap.get_arbitrary_string_data_from_pose(pose)`
        Data override precedences as defined in `ScoreMap::add_arbitrary_score_data_from_pose`:
            1. Per-residue real metrics
            2. Composite real metrics
            3. Real metrics
        Data override precedences as defined in `ScoreMap::add_arbitrary_string_data_from_pose`:
            1. Per-residue probabilities metrics
            2. Per-residue string metrics
            3. Composite string metrics
            4. String metrics
        """
        sm_string = self.string # Lowest precedence
        sm_composite_string = self.composite_string
        sm_per_residue_string = self.per_residue_string
        sm_per_residue_probabilities = self.per_residue_probabilities
        sm_real = self.real
        sm_composite_real = self.composite_real
        sm_per_residue_real = self.per_residue_real # Highest precedence

        _msg = "SimpleMetric {0} data key is clobbering SimpleMetric {1} data key: '{2}'"
        for k in sm_string.keys():
            if k in sm_composite_string.keys():
                self._clobber_warning(_msg.format("composite string", "string", k))
            if k in sm_per_residue_string.keys():
                self._clobber_warning(_msg.format("per-residue string", "string", k))
            if k in sm_per_residue_probabilities.keys():
                self._clobber_warning(_msg.format("per-residue probabilities", "string", k))
            if k in sm_real.keys():
                self._clobber_warning(_msg.format("real", "string", k))
            if k in sm_composite_real.keys():
                self._clobber_warning(_msg.format("composite real", "string", k))
            if k in sm_per_residue_real.keys():
                self._clobber_warning(_msg.format("per-residue real", "string", k))
        for k in sm_composite_string.keys():
            if k in sm_per_residue_string.keys():
                self._clobber_warning(_msg.format("per-residue string", "composite string", k))
            if k in sm_per_residue_probabilities.keys():
                self._clobber_warning(_msg.format("per-residue probabilities", "composite string", k))
            if k in sm_real.keys():
                self._clobber_warning(_msg.format("real", "composite string", k))
            if k in sm_composite_real.keys():
                self._clobber_warning(_msg.format("composite real", "composite string", k))
            if k in sm_per_residue_real.keys():
                self._clobber_warning(_msg.format("per-residue real", "composite string", k))
        for k in sm_per_residue_string.keys():
            if k in sm_per_residue_probabilities.keys():
                self._clobber_warning(_msg.format("per-residue probabilities", "per-residue string", k))
            if k in sm_real.keys():
                self._clobber_warning(_msg.format("real", "per-residue string", k))
            if k in sm_composite_real.keys():
                self._clobber_warning(_msg.format("composite real", "per-residue string", k))
            if k in sm_per_residue_real.keys():
                self._clobber_warning(_msg.format("per-residue real", "per-residue string", k))
        for k in sm_per_residue_probabilities.keys():
            if k in sm_real.keys():
                self._clobber_warning(_msg.format("real", "per-residue probabilities", k))
            if k in sm_composite_real.keys():
                self._clobber_warning(_msg.format("composite real", "per-residue probabilities", k))
            if k in sm_per_residue_real.keys():
                self._clobber_warning(_msg.format("per-residue real", "per-residue probabilities", k))
        for k in sm_real.keys():
            if k in sm_composite_real.keys():
                self._clobber_warning(_msg.format("composite real", "real", k))
            if k in sm_per_residue_real.keys():
                self._clobber_warning(_msg.format("per-residue real", "real", k))
        for k in sm_composite_real.keys():
            if k in sm_per_residue_real.keys():
                self._clobber_warning(_msg.format("per-residue real", "composite real", k))

        return types.MappingProxyType(
            {
                **sm_string,
                **sm_composite_string,
                **sm_per_residue_string,
                **sm_per_residue_probabilities,
                **sm_real,
                **sm_composite_real,
                **sm_per_residue_real,
            }
        )

    @property
    def string(self):
        return SimpleMetricStringDataAccessor(self.pose)

    @property
    def real(self):
        return SimpleMetricRealDataAccessor(self.pose)

    @property
    def composite_string(self):
        return SimpleMetricCompositeStringDataAccessor(self.pose)

    @property
    def composite_real(self):
        return SimpleMetricCompositeRealDataAccessor(self.pose)

    @property
    def per_residue_string(self):
        return SimpleMetricPerResidueStringDataAccessor(self.pose)

    @property
    def per_residue_real(self):
        return SimpleMetricPerResidueRealDataAccessor(self.pose)

    @property
    def per_residue_probabilities(self):
        return SimpleMetricPerResidueProbabilitiesDataAccessor(self.pose)

    def clear(self):
        """Clear pose SimpleMetric data."""
        clear_sm_data(self.pose)
