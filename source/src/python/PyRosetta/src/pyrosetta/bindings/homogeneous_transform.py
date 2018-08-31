from pyrosetta.bindings.utility import bind_method, bind_classmethod
from pyrosetta.rosetta.numeric import HomogeneousTransform_Double as HomogeneousTransform


# function taken from Will Sheffler's homog repo
def is_homog_xform(xforms):
    return ((xforms.shape[-2:] == (4, 4))
            and (np.allclose(1, np.linalg.det(xforms[..., :3, :3])))
            and (np.allclose(xforms[..., 3, :], [0, 0, 0, 1])))


@bind_classmethod(HomogeneousTransform)
def from_RT(cls, RTinstance):
    from pyrosetta.rosetta.core.kinematics import RT
    if not isinstance(RTinstance, RT):
        return
    return cls(RTinstance.get_rotation(), RTinstance.get_translation())


@bind_classmethod(HomogeneousTransform)
def from_array(cls, xform):
    from pyrosetta.rosetta.numeric import xyzVector_double_t as xyzVector
    import numpy as np
    if not isinstance(xform, np.ndarray) or not is_homog_xform(xform):
        return

    rot = xform[:3, :3]
    trans = xform[:3, 3]

    cols = [xyzVector(*rot[:, i]) for i in range(3)]
    cols.append(xyzVector(*trans))
    return cls(*cols)


@bind_method(HomogeneousTransform)
def as_array(self):
    """Convert self a numpy.array and return it.

    Returns:
        numpy.array: the 4x4 matrix representing the homogenous transform.
    """
    import numpy as np

    ht = np.identity(4)
    ht[:3, :3] = list(self.rotation_matrix())
    ht[:3, 3] = list(self.point())
    assert(is_homog_xform(ht))
    return ht
