# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# Pose bindings

import warnings

from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.pose import remove_upper_terminus_type_from_pose_residue
from pyrosetta.rosetta.core.pose import add_upper_terminus_type_to_pose_residue
from pyrosetta.rosetta.core.conformation import Residue
from pyrosetta.bindings.utility import slice_1base_indicies, bind_method, bind_property


class PoseResidueAccessor(object):
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
