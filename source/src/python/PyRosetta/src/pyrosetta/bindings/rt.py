from pyrosetta.bindings.utility import bind_method
from pyrosetta.rosetta.core.kinematics import RT


@bind_method(RT)
def as_array(self):
    """Convert self a numpy.array and return it.

    Returns:
        numpy.array: the 4x4 matrix representing the homogenous transform.
    """
    import numpy as np

    ht = np.identity(4)
    ht[:3, :3] = list(self.get_rotation())
    ht[:3, 3] = list(self.get_translation())
    assert(is_homog_xform(ht))
    return ht