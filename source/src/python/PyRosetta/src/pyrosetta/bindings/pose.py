# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# Pose bindings

import base64
import collections
import itertools
import pickle
import warnings


try:
    from collections.abc import Sized, Iterable, MutableSet, MutableMapping, Mapping
except ImportError:
    # these types are in the collections module in python < 3.3
    from collections import Sized, Iterable, MutableSet, MutableMapping, Mapping

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

from pyrosetta.rosetta.core.io.raw_data import ScoreMap
from pyrosetta.rosetta.core.simple_metrics import clear_sm_data, get_sm_data
from pyrosetta.rosetta.core.simple_metrics.metrics import (
    CustomRealValueMetric,
    CustomStringValueMetric,
)
from pyrosetta.rosetta.core.pose import remove_upper_terminus_type_from_pose_residue
from pyrosetta.rosetta.core.pose import add_upper_terminus_type_to_pose_residue
from pyrosetta.rosetta.core.conformation import Residue
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.bindings.utility import slice_1base_indicies, bind_method, bind_property
from pyrosetta.distributed.utility.pickle import (
    __cereal_getstate__,
    __cereal_setstate__,
)

from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.select.residue_selector import (
    TrueResidueSelector,
    NeighborhoodResidueSelector,
)
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueIndexSelector,
    AndResidueSelector,
)
from pyrosetta.rosetta.core.select.residue_selector import (
    PrimarySequenceNeighborhoodSelector,
)
from pyrosetta.rosetta.core.select.residue_selector import NotResidueSelector

def __pose_getstate__(pose):
    if pose.constraint_set().n_sequence_constraints() > 0:
        warnings.warn(
            "Pose pickling does not support sequence constraints.",
            UserWarning,
            stacklevel=2,
        )

        work_pose = Pose()
        work_pose.detached_copy(pose)
        work_pose.constraint_set().clear_sequence_constraints()

        return __pose_getstate__(work_pose)

    return __cereal_getstate__(pose)


Pose.__getstate__ = __pose_getstate__
Pose.__setstate__ = __cereal_setstate__


@bind_method(Pose)
def residue_pairs(
    self,
    primary_residue_selector=TrueResidueSelector(),
    secondary_residue_selector=TrueResidueSelector(),
    neighborhood_distance_maximum=None,
    sequence_distance_minimum=None,
):
    """Iterate the permutations of residue pairs from two selections which are
    within a cartesian and/or outside a sequence distance cutoff.

    The method uses the PrimarySequenceNeighborhoodSelector and
    NeighborhoodResidueSelector under the hood. Note that all _permutations_ of
    residues are yielded, not _combinations_. If you prefer combinations simply
    check `if residue_pair[0].sequence_position() > residue_pair[1].sequence_position():`.

    primary_residue_selector - ResidueSelector
    secondary_residue_selector - ResidueSelector
    neighborhood_distance - float
    sequence_distance - int

    return - iterator(tuple(Residue, Residue))
    """
    primary_residue_selection = primary_residue_selector.apply(self)
    primary_residue_indices = get_residues_from_subset(primary_residue_selection)

    for primary_residue_index in primary_residue_indices:
        temp_secondary_residue_selector = secondary_residue_selector

        primary_residue_index_selector = ResidueIndexSelector(primary_residue_index)
        temp_secondary_residue_selector = AndResidueSelector(
            temp_secondary_residue_selector,
            NotResidueSelector(primary_residue_index_selector),
        )

        if (
            neighborhood_distance_maximum
        ):  # Select residues within cartesian neighborhood distance
            neighborhood_residue_selector = NeighborhoodResidueSelector(
                primary_residue_index_selector, neighborhood_distance_maximum, False
            )
            temp_secondary_residue_selector = AndResidueSelector(
                temp_secondary_residue_selector, neighborhood_residue_selector
            )

        if (
            sequence_distance_minimum
        ):  # Select residues outside sequence neighborhood distance
            sequence_residue_selector = PrimarySequenceNeighborhoodSelector(
                sequence_distance_minimum,
                sequence_distance_minimum,
                primary_residue_index_selector,
            )
            temp_secondary_residue_selector = AndResidueSelector(
                temp_secondary_residue_selector,
                NotResidueSelector(sequence_residue_selector),
            )

        secondary_residue_selection = temp_secondary_residue_selector.apply(self)
        secondary_residue_indices = get_residues_from_subset(
            secondary_residue_selection
        )

        for secondary_residue_index in secondary_residue_indices:
            yield tuple(
                [
                    self.residues[primary_residue_index],
                    self.residues[secondary_residue_index],
                ]
            )


class PoseResidueAccessor(Sized, Iterable):
    """Accessor wrapper for pose objects providing a collection-like interface."""

    def __init__(self, pose):
        self.pose = pose

    def __len__(self):
        """The number of residues in pose."""
        return self.pose.size()

    def __iter__(self):
        """Iterate over the residues within pose."""
        for i in range(1, len(self) + 1):
            yield self.pose.residue(i)

    def __reversed__(self):
        """Iterate over the residues within pose in reverse order.
        Must provide __reversed__ as accessor object does not support 0-based
        indexing for the sequence protocol.
        """

        for i in range(len(self), 0, -1):
            yield self.pose.residue(i)

    def __getitem__(self, key):
        """1-based index and slice over residues.

        Examples: pose.residues[3], pose.residues[2:6], pose.residues[:-3],
                  pose.residues[2:], pose.residues[-1]
        """
        if isinstance(key, slice):
            return (self[i] for i in range(*slice_1base_indicies(key, len(self))))
        else:
            if key == 0:
                raise IndexError("1 base indexing does not support 0 index")
            if key < 0:
                key = len(self) + 1 + key
            return self.pose.residue(key)

    def __iadd__(self, other):
        """
        Short notation for appending a residue to the end of a pose by bond. If the
        terminal residue has Cterm variant type, the new appended residue will be
        transferred the variant type. If the pose is empty, the residue will be
        appended by jump for ease of use.
        """
        if isinstance(other, Residue):
            if not other.is_polymer():
                raise ValueError(
                    "Appending by bond on non-polymer, use *= for appending by jump"
                )

            if not len(self):
                self.pose.append_residue_by_jump(other, 1)
            else:
                if other.is_lower_terminus():
                    raise ValueError(
                        "Residue must not have lower terminus variant type."
                    )

                if self[-1].is_upper_terminus():
                    remove_upper_terminus_type_from_pose_residue(self.pose, len(self))
                    self.pose.append_polymer_residue_after_seqpos(
                        other, len(self), True
                    )
                    add_upper_terminus_type_to_pose_residue(self.pose, len(self))
                else:
                    self.pose.append_polymer_residue_after_seqpos(
                        other, len(self), True
                    )
        else:
            raise ValueError(
                "Only Residue objects can be appended to Pose residue accessor."
            )

        return self

    def __imul__(self, other):
        """
        Short notation for appending residues by jump to the terminal residue of a
        pose.
        """
        if isinstance(other, Residue):
            if not len(self):
                self.pose.append_residue_by_jump(other, 1)
            else:
                self.pose.append_residue_by_jump(other, len(self))
        else:
            raise ValueError(
                "Only Residue objects can be appended to Pose residue accessor."
            )
        return self


@bind_property(Pose)
def __residue_accessor(self):
    return PoseResidueAccessor(self)


@Pose.__residue_accessor.setter
def __residue_accessor(self, accessor):
    if not isinstance(accessor, PoseResidueAccessor) or accessor.pose is not self:
        raise AttributeError("Can't set residue accessor.")
    pass


Pose.residues = __residue_accessor


class ResidueLabelAccessor(MutableSet):
    """Accessor wrapper for a single residue's labels providing a mutable-set interface."""

    __slots__ = ("pdb_info", "res")

    def __init__(self, pdb_info, res):
        self.pdb_info = pdb_info
        self.res = res

    def __contains__(self, label):
        return self.pdb_info.res_haslabel(self.res, label)

    def __iter__(self):
        return iter(self.pdb_info.get_reslabels(self.res))

    def __len__(self):
        return len(self.pdb_info.get_reslabels(self.res))

    def add(self, label):
        self.pdb_info.add_reslabel(self.res, label)

    def clear(self):
        self.pdb_info.clear_reslabel(self.res)

    def discard(self, label):
        labels = set(self)
        if label not in labels:
            return

        self.clear()
        labels.discard(label)
        for l in labels:
            self.add(l)

    def __str__(self):
        return str(set(self))

    def _repr_pretty_(self, p, cycle):
        """IPython-display representation."""
        p.pretty(set(self))

    def __repr__(self):
        return "ResidueLabelAccessor(pdb_info=%r, res=%r)" % (self.pdb_info, self.res)


class PoseResidueLabelMaskAccessor(Mapping):
    """Accessor wrapper for residue label masks."""

    __slots__ = ("label_accessor",)

    def __init__(self, label_accessor):
        self.label_accessor = label_accessor

    def keys(self):
        return frozenset(itertools.chain.from_iterable(self.label_accessor))

    def __getitem__(self, label):
        return [label in r for r in self.label_accessor]

    def __len__(self):
        return len(self.keys())

    def __iter__(self):
        return iter(self.keys())


class PoseResidueLabelAccessor(object):
    """Accessor wrapper for a pose's residue labels as a 1-indexed-collection of sets."""

    __slots__ = ("pose",)

    @property
    def mask(self):
        return PoseResidueLabelMaskAccessor(self)

    def __init__(self, pose):
        self.pose = pose

    def __len__(self):
        return self.pose.pdb_info().nres()

    def __iter__(self):
        return (self[res] for res in range(1, len(self) + 1))

    def __reversed__(self):
        return (self[res] for res in range(len(self), 0, -1))

    def __getitem__(self, key):
        """1-based index and slice over residue labels."""
        if isinstance(key, slice):
            return (self[i] for i in range(*slice_1base_indicies(key, len(self))))
        else:
            if key == 0:
                raise IndexError("1 base indexing does not support 0 index")
            if key < 0:
                key = len(self) + 1 + key
            return ResidueLabelAccessor(self.pose.pdb_info(), key)

    def __str__(self):
        return list(map(set, self))

    def _repr_pretty_(self, p, cycle):
        """IPython-display representation."""
        p.pretty(list(self))

    def __repr__(self):
        return "PoseResidueLabelAccessor(pose=%s)" % self.pose


@bind_property(Pose)  # noqa: F811
def __reslabels_accessor(self):
    return PoseResidueLabelAccessor(self)


@Pose.__reslabels_accessor.setter
def __reslabels_accessor(self, accessor):
    if not isinstance(accessor, PoseResidueLabelAccessor) or accessor.pose is not self:
        raise AttributeError("Can't set reslabels accessor.")
    pass


Pose.reslabels = __reslabels_accessor


class PoseScoreSerializerBase(object):
    """Base class for `PoseScoreSerializer` methods."""
    __author__ = "Jason C. Klima"

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
        # Automatically serialize an arbitrary score value:
        pose.scores["foo"] = value

        # Manually serialize an arbitrary score value:
        value = PoseScoreSerializer.maybe_encode(value)

        # Automatically deserialize an arbitrary score value:
        value = pose.scores["foo"]

        # Manually deserialize an arbitrary score value:
        value = PoseScoreSerializer.maybe_decode(value)
    """
    __author__ = "Jason C. Klima"

    # Define different data types with human-readable custom prefixes in case anyone
    # accesses the serialized score values outside the scope of the `PoseScoreSerializer`
    _CustomTypeMetric = collections.namedtuple(
        "CustomTypeMetric", ["type", "prefix", "encode_func", "decode_func"],
    )
    _reserved_types = (str, float)
    _custom_type_metrics = {
        "bool": _CustomTypeMetric(
            type=bool,
            prefix="[CustomBoolMetric]",
            encode_func=str,
            decode_func=PoseScoreSerializerBase.bool_from_str,
        ),
        "int": _CustomTypeMetric(
            type=int,
            prefix="[CustomIntMetric]",
            encode_func=str,
            decode_func=int,
        ),
        "bytes": _CustomTypeMetric(
            type=bytes,
            prefix="[CustomBytesMetric]",
            encode_func=PoseScoreSerializerBase.to_base64,
            decode_func=PoseScoreSerializerBase.from_base64,
        ),
        "object": _CustomTypeMetric(
            type=object,
            prefix="[CustomArbitraryMetric]",
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


class PoseScoreAccessor(MutableMapping, PoseScoreSerializer):
    """Accessor wrapper for pose energies and extra scores."""

    __slots__ = ("pose",)

    _reserved = {k for k in ScoreType.__dict__.keys() if not k.startswith("_")}


    def __init__(self, pose):
        self.pose = pose
        self.custom_real_value_metric = CustomRealValueMetric()
        self.custom_string_value_metric = CustomStringValueMetric()

    @property
    def extra(self):
        import types

        _arbitrary_string_data = dict(ScoreMap.get_arbitrary_string_data_from_pose(self.pose).items())
        _arbitrary_score_data = dict(ScoreMap.get_arbitrary_score_data_from_pose(self.pose).items())

        for k in _arbitrary_string_data.keys():
            if k in _arbitrary_score_data.keys():
                warnings.warn(
                    "Arbitrary score data key is clobbering arbitrary string data key: '{0}'".format(k)
                )

        return types.MappingProxyType(
            dict(
                list(_arbitrary_string_data.items())
                + list(_arbitrary_score_data.items())
            )
        )

    @property
    def energies(self):
        import types

        return types.MappingProxyType(self.pose.energies().active_total_energies())

    @property
    def all(self):
        import types

        _arbitrary_string_data = dict(ScoreMap.get_arbitrary_string_data_from_pose(self.pose).items())
        _arbitrary_score_data = dict(ScoreMap.get_arbitrary_score_data_from_pose(self.pose).items())
        _active_total_energies = dict(self.pose.energies().active_total_energies().items())

        for k in _arbitrary_string_data.keys():
            if k in _arbitrary_score_data.keys():
                warnings.warn(
                    "Arbitrary score data key is clobbering arbitrary string data key: '{0}'".format(k)
                )
            if k in _active_total_energies.keys():
                warnings.warn(
                    "Active total energy score key is clobbering arbitrary string data key: '{0}'".format(k)
                )
        for k in _arbitrary_score_data.keys():
            if k in _active_total_energies.keys():
                warnings.warn(
                    "Active total energy score key is clobbering arbitrary score data key: '{0}'".format(k)
                )

        return types.MappingProxyType(
            dict(
                list(_arbitrary_string_data.items())
                + list(_arbitrary_score_data.items())
                + list(_active_total_energies.items())
            )
        )

    @property
    def _sm_data_attrs(self):
        """Supported `SimpleMetricStruct` attributes to reset after `clear_sm_data`."""
        return (
            'real_data_',
            'string_data_',
        )

    @property
    def _unsupported_sm_data_attrs(self):
        """Unsupported `SimpleMetricStruct` attributes that cannot be reset after `clear_sm_data`."""
        return (
            'composite_real_data_',
            'composite_string_data_',
            'per_residue_real_data_',
            'per_residue_real_output_',
            'per_residue_string_data_',
            'per_residue_string_output_',
        )

    def __len__(self):
        return len(self.all)

    def __iter__(self):
        return iter(self.all)

    def __getitem__(self, key):
        return self.maybe_decode(self.all[key])

    def __setitem__(self, key, value):
        if key in self._reserved:
            raise ValueError(
                "Can not set score key with reserved energy name: %r" % key
            )

        if isinstance(value, float):
            m = self.custom_real_value_metric
        else:
            m = self.custom_string_value_metric
        m.set_value(self.maybe_encode(value))
        m.apply(out_label=key, pose=self.pose, override_existing_data=True)

    def __delitem__(self, key):
        if key in self._reserved:
            raise ValueError(
                "Can not delete score key with reserved energy name: %r "
                "Consider 'pose.scores.clear()' or 'pose.energies().clear()'" % key
            )

        sm_data = get_sm_data(self.pose).get_all_sm_data()
        if hasPoseExtraScore(self.pose, key) or hasPoseExtraScore_str(self.pose, key):
            clearPoseExtraScore(self.pose, key)
        elif sm_data.has_data():
            # SimpleMetric data cannot be mutated by anything other than a SimpleMetric.
            self._maybe_delete_key_from_sm_data(key, sm_data)
        else:
            raise KeyError(key)

    def clear(self):
        """ Clear pose energies, extra scores, and SimpleMetric data"""
        self.pose.energies().clear()
        clearPoseExtraScores(self.pose)
        clear_sm_data(self.pose)

    def _get_sm_data(self, sm_data, attributes):
        """Get `SimpleMetricStruct` attributes as a dictionary."""
        return {
            _attr: dict(getattr(sm_data, _attr))
            for _attr in attributes
        }

    def _maybe_delete_key_from_sm_data(self, key, sm_data):
        """
        Cache, clear, and restore all SimpleMetric data except the deleted key.
        This is necessary to delete a single key, since `clear_sm_data` clears all
        SimpleMetric data.
        """
        # Cache SimpleMetric data
        _sm_data_dict = self._get_sm_data(sm_data, self._sm_data_attrs)
        _unsupported_sm_data_dict = self._get_sm_data(sm_data, self._unsupported_sm_data_attrs)
        # Raise if deleting an item we cannot restore
        for _attr, _d in _unsupported_sm_data_dict.items():
            for _k, _v in _d.items():
                if key == _k:
                    raise KeyError(
                        "'{0}'. Cannot delete a SimpleMetric data item from: {1}. ".format(
                            key, self._unsupported_sm_data_attrs
                        ) + "Consider using `pose.scores.clear()` or "
                        + "`pyrosetta.rosetta.core.simple_metrics.clear_sm_data(pose)`."
                    )
        if any(map(len, _unsupported_sm_data_dict.values())):
            raise KeyError(
                "'{0}'. Cannot delete item because Pose contains SimpleMetric data from: {1}".format(
                    key, self._unsupported_sm_data_attrs
                ) + "Consider using `pose.scores.clear()` or "
                + "`pyrosetta.rosetta.core.simple_metrics.clear_sm_data(pose)`."
            )
        # Clear SimpleMetric data
        clear_sm_data(self.pose)
        # Reapply all SimpleMetric data except the deleted item
        for _attr, _d in _sm_data_dict.items():
            for _k, _v in _d.items():
                if _k != key:
                    if _attr == "real_data_":
                        m = self.custom_real_value_metric
                    elif _attr == "string_data_":
                        m = self.custom_string_value_metric
                    else:
                        raise NotImplementedError(_attr)
                    m.set_value(_v)
                    m.apply(out_label=_k, pose=self.pose, override_existing_data=True)

    def __str__(self):
        return str(dict(self))

    def _repr_pretty_(self, p, cycle):
        """IPython-display representation."""
        p.pretty(dict(self))


@bind_property(Pose)  # noqa: F811
def __scores_accessor(self):
    return PoseScoreAccessor(self)


@Pose.__scores_accessor.setter
def __scores_accessor(self, accessor):
    if not isinstance(accessor, PoseScoreAccessor) or accessor.pose is not self:
        raise AttributeError("Can't set scores accessor.")
    pass


Pose.scores = __scores_accessor

# Deprecated Bindings - 19/11/17
@bind_method(Pose)
def __getitem__(self, key):
    warnings.warn(
        "Pose.__getitem__ is deprecated, prefer 'pose.residues.__getitem__'.",
        DeprecationWarning,
        stacklevel=2,
    )
    return self.residues.__getitem__(key)


@bind_method(Pose)
def __iter__(self):
    warnings.warn(
        "Pose.__iter__ is deprecated, prefer 'pose.residues.__iter__'.",
        DeprecationWarning,
        stacklevel=2,
    )
    return self.residues.__iter__()


@bind_method(Pose)
def __len__(self):
    warnings.warn(
        "Pose.__len__ is deprecated, prefer 'pose.residues.__len__'.",
        DeprecationWarning,
        stacklevel=2,
    )
    return len(self.residues)


@bind_method(Pose)
def pdb_rsd(self, chain_and_resNo):
    """Look up a specific PDB-numbered residue and return it.

    Args:
        chain_and_resNo (tuple): a tuple representing the PDB description of the residue
            in the format (chainID, resNo). For example, residue 1 on chain A would be
            ("A", 1).

    Returns:
        pyrosetta.core.conformation.Residue or None: the Residue instance in the Pose.
        returns `None` if the PDB residue identifier is invalid.
    """
    try:
        return self.residues[self.pdb_info().pdb2pose(*chain_and_resNo)]
    except IndexError:
        return


@bind_method(Pose)
def apply_transform(self, xform):
    """Apply a homogeneous transform to the current pose.

    Args:
        xform (np.ndarray): A homogeneous transform.
    """
    from pyrosetta.bindings.homogeneous_transform import is_homog_xform
    from pyrosetta.rosetta.numeric import (
        HomogeneousTransform_Double as HomogeneousTransform,
    )

    assert is_homog_xform(xform)

    ht = HomogeneousTransform.from_array(xform)
    self.apply_transform_Rx_plus_v(ht.rotation_matrix(), ht.point())


@bind_method(Pose)
def translate(self, t):
    """Apply a translation to all of the coordinates in a Pose.

    Args:
        p (Pose): The Pose instance to manipulate
        t (np.array): A vector to add to the Pose coordinates
    """
    import numpy as np

    assert t.shape in ((3,), (4,))
    xform = np.identity(4)
    xform[: t.shape[0], 3] = t
    self.apply_transform(xform)


@bind_method(Pose)
def rotate(self, R):
    """Apply a rotation matrix to all of the coordinates in a Pose.

    Args:
        p (Pose): The Pose instance to manipulate
        R (np.mat): A rotation matrix to apply to the Pose coordinates
    """
    import numpy as np

    assert R.shape in ((3, 3), (4, 3))
    xform = np.identity(4)
    xform[: R.shape[0], :3] = R
    self.apply_transform(xform)


@bind_method(Pose)
def get_hbonds(
    self,
    calculate_derivative=False,
    exclude_bb=False,
    exclude_bsc=False,
    exclude_scb=False,
    exclude_sc=False,
):
    """Return the HBondSet for all hydrogen bonds in the Pose.

    Args:
        caclulate_derivative (bool): Defaults to False.
        exclube_bb (bool): Exclude backbone hydrogen bonds from the returned HBondSet.
            Defaults to False.
        exclube_bsc (bool):  Exclude backbone--side-chain hydrogen bonds from the returned HBondSet.
            Defaults to False.
        exclube_scb (bool):  Exclude side-chain--backbone hydrogen bonds from the returned HBondSet.
            Defaults to False.
        exclube_sc (bool):  Exclude side-chain hydrogen bonds from the returned HBondSet.
            Defaults to False.

    Returns:
        pyrosetta.rosetta.core.scoring.hbonds.HBondSet: THe selected hydrogen bonds in the Pose.
    """
    self.update_residue_neighbors()
    import pyrosetta.rosetta.core.scoring.hbonds

    hbset = pyrosetta.rosetta.core.scoring.hbonds.HBondSet()
    pyrosetta.rosetta.core.scoring.hbonds.fill_hbond_set(
        self,
        calculate_derivative,
        hbset,
        exclude_bb,
        exclude_bsc,
        exclude_scb,
        exclude_sc,
    )
    return hbset


@bind_method(Pose)
def display_secstruct(self, space=8, page=80):
    """Use DSSP to assign secondary structure codes to each residue in a Pose,
    store the result in the Pose, and display the result.

    Args:
        space (int): Spacing of residue number labels in display.
            Defaults to 8.
        page (int): Number of characters per line in dispaly.
            Defaults to 80.
    """
    import pyrosetta.rosetta.core.scoring.dssp

    dssp = pyrosetta.rosetta.core.scoring.dssp.Dssp(self)
    dssp.insert_ss_into_pose(self)

    seq = self.sequence()
    sec = self.secstruct()
    assert len(seq) == len(sec)

    seq_struc = [
        (seq[i : i + page], sec[i : i + page]) for i in range(0, len(seq), page)
    ]
    for i, (seq, struc) in enumerate(seq_struc):
        num = "".join(
            [
                str(i).ljust(space)
                for i in range(page * i + 1, page * i + min(page, len(seq)), space)
            ]
        )
        print("\n".join([num, seq, struc, "\n"]))
