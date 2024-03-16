# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"


from pyrosetta.bindings.utility import bind_method
from pyrosetta.rosetta.core.select import get_residues_from_subset
from pyrosetta.rosetta.core.select.residue_selector import (
    AndResidueSelector,
    NotResidueSelector,
    OrResidueSelector,
    ResidueSelector,
)


@bind_method(ResidueSelector)
def get_residues(self, pose):
    """
    Return a python `list` object of selected residues in pose numbering.
    
    Args:
        pose: a `Pose` object to which to apply this residue selector.

    Returns:
        A `list` object of selected residues in pose numbering.
    """
    return list(get_residues_from_subset(self.apply(pose)))


@bind_method(ResidueSelector)
def __and__(self, selector):
    """
    Short notation for `AndResidueSelector` using the `&` operator in python.
    """
    return AndResidueSelector(self, selector)


@bind_method(ResidueSelector)
def __or__(self, selector):
    """
    Short notation for `OrResidueSelector` using the `|` operator in python.
    """
    return OrResidueSelector(self, selector)


@bind_method(ResidueSelector)
def __invert__(self):
    """
    Short notation for `NotResidueSelector` using the `~` operator in python.
    """
    return NotResidueSelector(self)


@bind_method(ResidueSelector)
def __xor__(self, selector):
    """
    Short notation for exclusive `OrResidueSelector` using the `^` operator in python.
    """
    return AndResidueSelector(OrResidueSelector(self, selector), NotResidueSelector(AndResidueSelector(self, selector)))


@bind_method(ResidueSelector)
def __iand__(self, selector):
    """
    Short notation for in-place `AndResidueSelector` using the `&=` operator in python.
    """
    return self.__and__(selector)


@bind_method(ResidueSelector)
def __ior__(self, selector):
    """
    Short notation for in-place `OrResidueSelector` using the `|=` operator in python.
    """
    return self.__or__(selector)


@bind_method(ResidueSelector)
def __ixor__(self, selector):
    """
    Short notation for in-place exclusive `OrResidueSelector` using the `^=` operator in python.
    """
    return self.__xor__(selector)

