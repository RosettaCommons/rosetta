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

from pyrosetta.bindings.scores.base import PoseCacheAccessorBase


class EnergiesAccessor(PoseCacheAccessorBase, MutableMapping):
    """Accessor wrapper for pose energies."""
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
