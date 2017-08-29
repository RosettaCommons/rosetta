# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# Pose dunder bindings
# @author atom-moyer

from pyrosetta.util import bind_method, vector1_indices
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.pose import pose_residue_is_terminal
from pyrosetta.rosetta.core.pose import remove_upper_terminus_type_from_pose_residue
from pyrosetta.rosetta.core.pose import remove_lower_terminus_type_from_pose_residue
from pyrosetta.rosetta.core.pose import add_upper_terminus_type_to_pose_residue
from pyrosetta.rosetta.core.conformation import Residue

@bind_method(Pose)
def __len__(self):
    """
    Define the dunder len function for pose as the amount of residues.

    @atom-moyer
    """
    return self.size()


@bind_method(Pose)
def __iter__(self):
    """
    Define the dunder iter function for pose as stream of ordered residues.

    @atom-moyer
    """
    for i in range(1, len(self)+1):
        yield self.residue(i)


@bind_method(Pose)
def __getitem__(self, key):
    """
    Quick index and slicing operation.  There is full support for
    conventional python style slicing.

    Examples: pose[3], pose[2:6], pose[:-3], pose[2:], pose[-1]

    @atom-moyer
    """
    if isinstance(key, slice): # Support for Vector1 Index Slicing
        key_vector1_indices = vector1_indices(key, len(self))
        return [self.residue(i) for i in range(*key_vector1_indices)]
    else:
        if key < 0: key = len(self) + 1 + key
        return self.residue(key)


@bind_method(Pose)
def __iadd__(self, other):
    """
    Short notation for appending a residue to the end of a pose by bond. If the
    terminal residue has Cterm variant type, the new appended residue will be
    transferred the variant type. If the pose is empty, the residue will be
    appended by jump for ease of use.

    @atom-moyer
    """
    if isinstance(other, Residue):
        assert other.is_polymer(), 'Only appending by bond is viable for polymer, use *= for appending by jump'
        if not len(self):
            print('Appending residue by jump because appended pose is empty.')
            self.append_residue_by_jump(other, 1)
        else:
            assert not other.is_lower_terminus(), 'Residue must not have lower terminus variant type.'
            if self[-1].is_upper_terminus():
                remove_upper_terminus_type_from_pose_residue(self, len(self))
                self.append_polymer_residue_after_seqpos(other, len(self), True)
                add_upper_terminus_type_to_pose_residue(self, len(self))
            else:
                self.append_polymer_residue_after_seqpos(other, len(self), True)
    else:
        raise ValueError('Only Residue objects can be appended to Pose objects')
    return self


@bind_method(Pose)
def __imul__(self, other):
    """
    Short notation for appending residues by jump to the terminal residue of a
    pose.

    @atom-moyer
    """
    if isinstance(other, Residue):
        if not len(self):
            self.append_residue_by_jump(other, 1)
        else:
            self.append_residue_by_jump(other, len(self))
    else:
        raise ValueError('Only Residue objects can be appended to Pose objects')
    return self
