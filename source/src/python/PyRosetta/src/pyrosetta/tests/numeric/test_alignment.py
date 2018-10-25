import unittest
import numpy

import pyrosetta.numeric.alignment.rmsd_calc as rmsd_calc


class TestAlignment(unittest.TestCase):

    @staticmethod
    def _rotation_matrix(axis, theta):
        from numpy import cross, eye
        from scipy.linalg import expm, norm
        return expm(cross(eye(3), axis/norm(axis)*theta))

    def test_superimpose(self):
        numpy.random.seed(1663)
        test_system = numpy.random.random((10, 3)) * 10
        test_rotation = self._rotation_matrix([1, 2, 3], 1.3)

        test_coords = numpy.empty((2, 10, 3))
        for i in range(10):
            test_coords[0, i] = test_system[i]
            test_coords[1, i] = numpy.dot(test_rotation, test_system[i])

        i_rmsd = rmsd_calc.coordinate_array_rmsd(test_coords[0], test_coords[1])
        numpy.testing.assert_allclose(i_rmsd, 0)

        aligned, s_rmsd = rmsd_calc.superimpose_coordinate_array(test_coords, test_coords[1], return_rmsd=True)
        numpy.testing.assert_allclose(s_rmsd, 0, atol=1e-6)
        for i in range(aligned.shape[0]):
            numpy.testing.assert_allclose(numpy.mean(aligned[i], axis=0), numpy.mean(test_coords[1], axis=0))
            numpy.testing.assert_allclose(aligned[i], test_coords[1])

    def test_rmsd_calc(self):
        test_coords = numpy.zeros((2, 2, 3))
        test_coords[0, 0, 0] = 1
        test_coords[0, 1, 0] = -1

        test_coords[1, 0, 1] = 2
        test_coords[1, 1, 1] = -2

        numpy.testing.assert_allclose(
            rmsd_calc.coordinate_array_broadcast_rmsd(test_coords, test_coords),
            numpy.array([[0, 1], [1, 0]])
        )

        numpy.testing.assert_allclose(
            rmsd_calc.coordinate_array_broadcast_rmsd(test_coords[0], test_coords),
            numpy.array([0, 1])
        )

        numpy.testing.assert_allclose(
            rmsd_calc.coordinate_array_broadcast_rmsd(test_coords, test_coords[0]),
            numpy.array([0, 1])
        )

        numpy.testing.assert_allclose(
            rmsd_calc.coordinate_array_broadcast_rmsd(test_coords[0], test_coords[1]),
            numpy.array([1])
        )

        numpy.testing.assert_allclose(
            rmsd_calc.coordinate_array_rmsd(test_coords, test_coords),
            numpy.array([0, 0])
        )

        numpy.testing.assert_allclose(
            rmsd_calc.coordinate_array_rmsd(test_coords[0], test_coords),
            numpy.array([0, 1])
        )

        numpy.testing.assert_allclose(
            rmsd_calc.coordinate_array_rmsd(test_coords[[1, 0]], test_coords),
            numpy.array([1, 1])
        )
