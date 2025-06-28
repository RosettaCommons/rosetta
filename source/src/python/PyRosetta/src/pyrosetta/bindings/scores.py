# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# Pose bindings for scores


__author__ = "Jason C. Klima"


import base64
import collections
import pickle
import types
import warnings

try:
    from collections.abc import MutableMapping
except ImportError:
    # For python < 3.3
    from collections import MutableMapping

from ..rosetta.core.pose import Pose
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
from pyrosetta.rosetta.core.simple_metrics.metrics import (
    CustomRealValueMetric,
    CustomStringValueMetric,
)
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.bindings.utility import bind_property


class PoseScoreSerializerBase(object):
    """Base class for `PoseScoreSerializer` methods."""

    @staticmethod
    def to_pickle(value):
        try:
            return pickle.dumps(value)
        except (TypeError, OverflowError, MemoryError, pickle.PicklingError) as ex:
            raise TypeError(
                "Only pickle-serializable object types are allowed to be set "
                + "as score values. Received: %r. %s" % (type(value), ex)
            )

    @staticmethod
    def from_pickle(value):
        try:
            return pickle.loads(value)
        except (TypeError, OverflowError, MemoryError, EOFError, pickle.UnpicklingError) as ex:
            raise TypeError(
                "Could not deserialize score value of type %r. %s" % (type(value), ex)
            )

    @staticmethod
    def to_base64(value):
        return base64.b64encode(value).decode()

    @staticmethod
    def from_base64(value):
        return base64.b64decode(value)

    @staticmethod
    def to_base64_pickle(value):
        return PoseScoreSerializerBase.to_base64(PoseScoreSerializerBase.to_pickle(value))

    @staticmethod
    def from_base64_pickle(value):
        return PoseScoreSerializerBase.from_pickle(PoseScoreSerializerBase.from_base64(value))

    @staticmethod
    def bool_from_str(value):
        if value == "True":
            return True
        elif value == "False":
            return False
        else:
            raise NotImplementedError(value)


class PoseScoreSerializer(PoseScoreSerializerBase):
    """
    Serialize and deserialize score values for CustomStringValueMetric SimpleMetric.

    Examples:
        Automatically serialize an arbitrary score value:
            `pose.cache["foo"] = value`
        Automatically deserialize an arbitrary score value:
            `value = pose.cache["foo"]`
        Manually serialize an arbitrary score value:
            `value = PoseScoreSerializer.maybe_encode(value)`
        Manually deserialize an arbitrary score value:
            `value = PoseScoreSerializer.maybe_decode(value)`
    """
    # Define different data types with human-readable custom prefixes in case anyone
    # accesses the serialized score values outside the scope of the `PoseScoreSerializer`
    _CustomTypeMetric = collections.namedtuple(
        "CustomTypeMetric", ["type", "prefix", "encode_func", "decode_func"],
    )
    _reserved_types = (str, float)
    _custom_type_metrics = {
        "bool": _CustomTypeMetric(
            type=bool,
            prefix="[CustomBooleanValueMetric]",
            encode_func=str,
            decode_func=PoseScoreSerializerBase.bool_from_str,
        ),
        "int": _CustomTypeMetric(
            type=int,
            prefix="[CustomDiscreteValueMetric]",
            encode_func=str,
            decode_func=int,
        ),
        "bytes": _CustomTypeMetric(
            type=bytes,
            prefix="[CustomBinaryValueMetric]",
            encode_func=PoseScoreSerializerBase.to_base64,
            decode_func=PoseScoreSerializerBase.from_base64,
        ),
        "object": _CustomTypeMetric(
            type=object,
            prefix="[CustomArbitraryValueMetric]",
            encode_func=PoseScoreSerializerBase.to_base64_pickle,
            decode_func=PoseScoreSerializerBase.from_base64_pickle,
        ),
    }

    @staticmethod
    def maybe_encode(value):
        """Serialize the input value into a `str` object if it's not a `str` or `float` object."""
        if not isinstance(value, PoseScoreSerializer._reserved_types):
            for _custom_type_metric in PoseScoreSerializer._custom_type_metrics.values():
                if isinstance(value, _custom_type_metric.type):
                    value = "{0}{1}".format(
                        _custom_type_metric.prefix,
                        _custom_type_metric.encode_func(value)
                    )
                    break
            else:
                raise NotImplementedError("Unsupported object type: {0}".format(type(value)))

        return value

    @staticmethod
    def maybe_decode(value):
        """Deserialize the input value if it's serialized."""
        if isinstance(value, str):
            for _custom_type_metric in PoseScoreSerializer._custom_type_metrics.values():
                if value.startswith(_custom_type_metric.prefix):
                    value = _custom_type_metric.decode_func(
                        value[len(_custom_type_metric.prefix):]
                    )
                    break
            else:
                pass # Return without decoding because string doesn't start with a prefix

        return value


class ClobberWarning(UserWarning):
    pass


class PoseCacheAccessorBase(PoseScoreSerializer):
    __slots__ = ("pose",)

    def __init__(self, pose):
        self.pose = pose
        self.custom_real_value_metric = CustomRealValueMetric()
        self.custom_string_value_metric = CustomStringValueMetric()

    @property
    def _default_custom_metric_keys(self):
        return {"custom_real_valued_metric", "custom_string_valued_metric"}

    @property
    def _reserved(self):
        _reserved = {k for k in ScoreType.__dict__.keys() if not k.startswith("_")}
        _reserved = _reserved.union(self._default_custom_metric_keys)
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
        _accessor = SimpleMetricDataAccessor(self.pose)
        return {
            _attr: dict(getattr(_accessor, _attr))
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


class SimpleMetricDataAccessorBase(PoseCacheAccessorBase, MutableMapping):
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


class ExtraScoresAccessorBase(PoseCacheAccessorBase, MutableMapping):
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


class EnergiesAccessor(PoseCacheAccessorBase, MutableMapping):
    __slots__ = ("pose",)

    def __init__(self, pose):
        super().__init__(pose)

    @property
    def all(self):
        return types.MappingProxyType(
            dict(self.pose.energies().active_total_energies().items())
        )

    def __getitem__(self, key):
        return self.all[key]

    def __setitem__(self, key, value):
        raise NotImplementedError(
            "Cannot set an item to pose energies. Consider scoring the pose with a score function."
        )

    def __delitem__(self, key):
        raise NotImplementedError(
            "Cannot delete an item from pose energies. Consider using `pose.energies.clear()`."
        )

    def clear(self):
        """Clear pose energies."""
        self.pose.energies().clear()


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
                "extra": {
                    "string": dict(self.extra.string),
                    "real": dict(self.extra.real),
                },
                "metrics": {
                    "string": dict(self.metrics.string),
                    "real": dict(self.metrics.real),
                    "composite_string": dict(self.metrics.composite_string),
                    "composite_real": dict(self.metrics.composite_real),
                    "per_residue_string": dict(self.metrics.per_residue_string),
                    "per_residue_real": dict(self.metrics.per_residue_real),
                    "per_residue_probabilities": dict(self.metrics.per_residue_probabilities),
                },
                "energies": dict(self.energies),
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
        """Assert that `pose.cache.all_keys` returns all unique keys."""
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
        return self.maybe_decode(self.all[key])

    def __setitem__(self, key, value):
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


@bind_property(Pose)  # noqa: F811
def __cache_accessor(self):
    return PoseCacheAccessor(self)


@Pose.__cache_accessor.setter
def __cache_accessor(self, accessor):
    if not isinstance(accessor, PoseCacheAccessor) or accessor.pose is not self:
        raise AttributeError("Can't set cache accessor.")
    pass


Pose.cache = __cache_accessor
