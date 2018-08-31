# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# Pose bindings

import warnings

import itertools
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
    clearPoseExtraScore,
    clearPoseExtraScores,
)
from pyrosetta.rosetta.core.pose import remove_upper_terminus_type_from_pose_residue
from pyrosetta.rosetta.core.pose import add_upper_terminus_type_to_pose_residue
from pyrosetta.rosetta.core.conformation import Residue
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.bindings.utility import slice_1base_indicies, bind_method, bind_property
from pyrosetta.distributed.utility.pickle import __cereal_getstate__, __cereal_setstate__

from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector, NeighborhoodResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, AndResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import PrimarySequenceNeighborhoodSelector
from pyrosetta.rosetta.core.select.residue_selector import NotResidueSelector


def __pose_getstate__(pose):
    if pose.constraint_set().n_sequence_constraints() > 0:
        warnings.warn(
            "Pose pickling does not support sequence constraints.",
            UserWarning,
            stacklevel=2
        )

        work_pose = Pose()
        work_pose.detached_copy(pose)
        work_pose.constraint_set().clear_sequence_constraints()

        return __pose_getstate__(work_pose)

    return __cereal_getstate__(pose)


Pose.__getstate__ = __pose_getstate__
Pose.__setstate__ = __cereal_setstate__


@bind_method(Pose)
def residue_pairs(self,
                  primary_residue_selector=TrueResidueSelector(),
                  secondary_residue_selector=TrueResidueSelector(),
                  neighborhood_distance_maximum=None,
                  sequence_distance_minimum=None):
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
        temp_secondary_residue_selector = AndResidueSelector(temp_secondary_residue_selector,
                                                             NotResidueSelector(primary_residue_index_selector))

        if neighborhood_distance_maximum: # Select residues within cartesian neighborhood distance
            neighborhood_residue_selector = NeighborhoodResidueSelector(primary_residue_index_selector,
                                                                        neighborhood_distance_maximum,
                                                                        False)
            temp_secondary_residue_selector = AndResidueSelector(temp_secondary_residue_selector,
                                                                 neighborhood_residue_selector)

        if sequence_distance_minimum: # Select residues outside sequence neighborhood distance
            sequence_residue_selector = PrimarySequenceNeighborhoodSelector(sequence_distance_minimum,
                                                                            sequence_distance_minimum,
                                                                            primary_residue_index_selector)
            temp_secondary_residue_selector = AndResidueSelector(temp_secondary_residue_selector,
                                                                 NotResidueSelector(sequence_residue_selector))

        secondary_residue_selection = temp_secondary_residue_selector.apply(self)
        secondary_residue_indices = get_residues_from_subset(secondary_residue_selection)

        for secondary_residue_index in secondary_residue_indices:
            yield tuple([self.residues[primary_residue_index], self.residues[secondary_residue_index]])


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
                    'Appending by bond on non-polymer, use *= for appending by jump')

            if not len(self):
                self.pose.append_residue_by_jump(other, 1)
            else:
                if other.is_lower_terminus():
                    raise ValueError(
                        'Residue must not have lower terminus variant type.')

                if self[-1].is_upper_terminus():
                    remove_upper_terminus_type_from_pose_residue(
                        self.pose, len(self))
                    self.pose.append_polymer_residue_after_seqpos(
                        other, len(self), True)
                    add_upper_terminus_type_to_pose_residue(
                        self.pose, len(self))
                else:
                    self.pose.append_polymer_residue_after_seqpos(
                        other, len(self), True)
        else:
            raise ValueError(
                'Only Residue objects can be appended to Pose residue accessor.')

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
                'Only Residue objects can be appended to Pose residue accessor.')
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
        return frozenset(
            itertools.chain.from_iterable(self.label_accessor))

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

class PoseScoreAccessor(MutableMapping):
    """Accessor wrapper for pose energies and extra scores."""
    __slots__ = ("pose",)

    _reserved = {
        k for k in ScoreType.__dict__.keys() if not k.startswith("_")
    }

    def __init__(self, pose):
        self.pose = pose

    @property
    def extra(self):
        import types
        return types.MappingProxyType(dict(
            list(getPoseExtraFloatScores(self.pose).items()) +
            list(getPoseExtraStringScores(self.pose).items())
        ))

    @property
    def energies(self):
        import types
        return types.MappingProxyType(
            self.pose.energies().active_total_energies()
        )

    @property
    def all(self):
        import types
        return types.MappingProxyType(dict(
            list(getPoseExtraFloatScores(self.pose).items()) +
            list(getPoseExtraStringScores(self.pose).items()) +
            list(self.pose.energies().active_total_energies().items())
        ))

    def __len__(self):
        return len(self.all)

    def __iter__(self):
        return iter(self.all)

    def __getitem__(self, key):
        return self.all[key]

    def __setitem__(self, key, value):
        if key in self._reserved:
            raise ValueError(
                "Can not set score key with reserved energy name: %r" % key)

        # Bit of a two-step to deal with potential duplicate keys in the
        # score maps. First check if a key, of either type, exists. If so
        # *try* to set the extra score, triggering type conversion checking
        # etc...
        #
        # If set is successful then clear the score cache and set again,
        # eliminating any potential duplicate keys.
        had_score = key in self

        setPoseExtraScore(self.pose, key, value)
        if had_score:
            self.__delitem__(key)
            setPoseExtraScore(self.pose, key, value)

    def __delitem__(self, key):
        if key in self._reserved:
            raise ValueError(
                "Can not delete score key with reserved energy name: %r "
                "Consider 'pose.scores.clear()' or 'pose.energies().clear()'"
                % key
            )

        if not hasPoseExtraScore(self.pose, key):
            raise KeyError(key)

        clearPoseExtraScore(self.pose, key)

    def clear(self):
        self.pose.energies().clear()
        clearPoseExtraScores(self.pose)

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
        stacklevel=2
    )
    return self.residues.__getitem__(key)


@bind_method(Pose)
def __iter__(self):
    warnings.warn(
        "Pose.__iter__ is deprecated, prefer 'pose.residues.__iter__'.",
        DeprecationWarning,
        stacklevel=2
    )
    return self.residues.__iter__()


@bind_method(Pose)
def __len__(self):
    warnings.warn(
        "Pose.__len__ is deprecated, prefer 'pose.residues.__len__'.",
        DeprecationWarning,
        stacklevel=2
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
    from pyrosetta.rosetta.numeric import HomogeneousTransform_Double as HomogeneousTransform
    assert(is_homog_xform(xform))

    ht = HomogeneousTransform.from_array(xform)
    self.apply_transform_Rx_plus_v(ht.rotation_matrix(), ht.point())
