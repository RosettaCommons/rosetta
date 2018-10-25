import numpy

from pyrosetta.utility.array import structured_array_to_basic, basic_array_to_structured

import pyrosetta.rosetta.numeric.alignment


def superimpose_coordinate_array(src, onto, superimpose=None, rmsd_out=None, return_rmsd=False):
    """Calculate and apply superposition transform between src and onto coordinates.

    src - array(shape=([a], n, 3), dtype=float)
    onto- array(shape=([a], n, 3), dtype=float)
    superimpose - array(shape=([a], m, 3), dtype=float)
    rmsd_out - array(shape=(a), dtype=float)
    return_rmsd - bool

    src and onto must contain dimension n, the number of entries per coordinate set.

    superimpose may contain a different number of coordinates per entry, and
    will be updated in-place with coordinates updated by src->onto transform.

    returns:
        superimposed_coordinates - array(shape=(a, m, 3), dtype=float)
        OR
        (superimposed_coordinates, pairwise_rmsd) if return_rmsd
    """
    if superimpose is None:
        superimpose = src.copy()

    flatten = src.ndim == 2 and onto.ndim == 2

    if src.ndim == 2:
        src = numpy.expand_dims(src, 0)

    if onto.ndim == 2:
        onto = numpy.expand_dims(onto, 0)

    if superimpose.ndim == 2:
        superimpose = numpy.expand_dims(superimpose, 0)

    if rmsd_out is None:
        rmsd_out = numpy.empty(src.shape[0])
    else:
        rmsd_out = rmsd_out.reshape(src.shape[0])

    pyrosetta.rosetta.numeric.alignment.coordinate_array_superimpose(src, onto, superimpose, rmsd_out)

    if flatten:
        rmsd_out = rmsd_out.reshape(())
        superimpose = superimpose.reshape((-1, 3))

    if return_rmsd:
        return superimpose, rmsd_out
    else:
        return superimpose


def superimpose_structured_array(src, onto, superimpose=None):
    """Calculate and apply superposition transform between src and onto coordinates.

    src         - array(shape=([a], n), dtype=[(name, float, 3)...])
    onto        - array(shape=([a|1], n), dtype=[(name, float, 3)...])
    superimpose -   array(shape=([a], n_2), dtype=[(name, float, 3)...])

    src and onto must contain dimension n and be of the same structured dtype.

    returns superimposed - array(shape=([a], n_2), dtype=[(name, float, 3)...])
    """

    assert src.shape[-1] == onto.shape[-1]
    assert src.dtype == onto.dtype

    if src.ndim == 1:
        ac = structured_array_to_basic(src).reshape((-1, 3))
    else:
        assert src.ndim == 2
        ac = structured_array_to_basic(src).reshape((src.shape[0], -1, 3))

    if onto.ndim == 1:
        bc = structured_array_to_basic(onto).reshape((-1, 3))
    else:
        assert onto.ndim == 2
        bc = structured_array_to_basic(onto).reshape((onto.shape[0], -1, 3))

    if superimpose is None:
        superimpose = src.copy()

    if superimpose.ndim == 1:
        sc = structured_array_to_basic(superimpose).reshape((-1, 3))
    else:
        assert superimpose.ndim == 2
        sc = structured_array_to_basic(
            superimpose).reshape((onto.shape[0], -1, 3))

    superimpose_coordinate_array(ac, bc, sc)
    sc = sc.reshape(sc.shape[:-2] + (-1, len(superimpose.dtype.names), 3))

    structured_result = basic_array_to_structured(
        sc,
        superimpose.dtype.names,
        superimpose.dtype[superimpose.dtype.names[0]])

    return numpy.squeeze(structured_result, axis=-1)


def coordinate_array_rmsd(coordinates_a, coordinates_b, out=None):
    """Calculate rmsd between entries in the given coordinate arrays.

    coordinates_a - array(shape=([a], n, 3), dtype=float)
    coordinates_b - array(shape=([a], n, 3), dtype=float)
    out - array(shape=(a), dtype=float)

    coordinates_a and coordinates_b must contain dimension n, the number of entries per coordinate set.

    returns pairwise_rmsd - array(shape=(n), dtype=float) or shape () array-scalar if both inputs are 2-dim
    """

    flatten = coordinates_a.ndim == 2 and coordinates_b.ndim == 2

    if coordinates_a.ndim == 2:
        coordinates_a = numpy.expand_dims(coordinates_a, 0)

    if coordinates_b.ndim == 2:
        coordinates_b = numpy.expand_dims(coordinates_b, 0)

    if out is None:
        out = numpy.empty(max(coordinates_a.shape[0], coordinates_b.shape[0]))
    else:
        out = out.reshape(max(coordinates_a.shape[0], coordinates_b.shape[0]))

    pyrosetta.rosetta.numeric.alignment.coordinate_array_rmsd(
        coordinates_a, coordinates_b, out)

    if flatten:
        return out.reshape(())
    else:
        return out


def coordinate_array_broadcast_rmsd(coordinates_a, coordinates_b, out=None):
    """Calculate rmsd between entries in the given coordinate arrays.

    coordinates_a - array(shape=([a], n, 3), dtype=float)
    coordinates_b - array(shape=([b], n, 3), dtype=float)
    out - array(shape=(a, b), dtype=float)

    coordinates_a and coordinates_b must contain dimension n, the number of entries per coordinate set.

    returns broadcast_rmsd - array(shape=(a,b), dtype=float)
        Flattened to 1 dimension if a or b is dimension 2.
    """

    flatten = False

    if coordinates_a.ndim == 2:
        flatten = True
        coordinates_a = numpy.expand_dims(coordinates_a, 0)

    if coordinates_b.ndim == 2:
        flatten = True

        coordinates_b = numpy.expand_dims(coordinates_b, 0)

    if out is None:
        out = numpy.empty((coordinates_a.shape[0], coordinates_b.shape[0]))
    else:
        out = out.reshape((coordinates_a.shape[0], coordinates_b.shape[0]))

    pyrosetta.rosetta.numeric.alignment.coordinate_array_broadcast_rmsd(
        coordinates_a, coordinates_b, out)

    if flatten:
        return out.ravel()
    else:
        return out


def structured_array_broadcast_rmsd(coordinates_a, coordinates_b):
    """Calculate rmsd between entries in the given structured coordinate arrays.

    coordinates_a - array(shape=([a], n), dtype=[(name, float, 3)...])
    coordinates_b - array(shape=([b], n), dtype=[(name, float, 3)...])
    out - array(shape=(a, b), dtype=float)

    coordinates_a and coordinates_b must contain dimension n and be of the same structured dtype.

    returns broadcast_rmsd - array(shape=(a,b), dtype=float)
        Flattened to 1 dimension if a or b is dimension 2.
    """
    assert coordinates_a.shape[-1] == coordinates_b.shape[-1]
    assert coordinates_a.dtype == coordinates_b.dtype

    if coordinates_a.ndim == 1:
        ac = structured_array_to_basic(coordinates_a).reshape((-1, 3))
    else:
        assert coordinates_a.ndim == 2
        ac = structured_array_to_basic(coordinates_a).reshape(
            (coordinates_a.shape[0], -1, 3))

    if coordinates_b.ndim == 1:
        bc = structured_array_to_basic(coordinates_b).reshape((-1, 3))
    else:
        assert coordinates_b.ndim == 2
        bc = structured_array_to_basic(coordinates_b).reshape(
            (coordinates_b.shape[0], -1, 3))

    return coordinate_array_broadcast_rmsd(ac, bc)


def structured_array_unaligned_broadcast_rmsd(coordinates_a, coordinates_b):
    """Calculate rmsd between entries in the given structured coordinate arrays.

    coordinates_a - array(shape=([a], n), dtype=[(name, float, 3)...])
    coordinates_b - array(shape=([b], n), dtype=[(name, float, 3)...])
    out - array(shape=(a, b), dtype=float)

    coordinates_a and coordinates_b must contain dimension n and be of the same structured dtype.

    returns broadcast_rmsd - array(shape=(a,b), dtype=float)
        Flattened to 1 dimension if a or b is dimension 2.
    """
    assert coordinates_a.shape[-1] == coordinates_b.shape[-1]
    assert coordinates_a.dtype == coordinates_b.dtype

    if coordinates_a.ndim == 1:
        ac = structured_array_to_basic(coordinates_a).reshape((-1, 3))
    else:
        assert coordinates_a.ndim == 2
        ac = structured_array_to_basic(coordinates_a).reshape(
            (coordinates_a.shape[0], -1, 3))

    if coordinates_b.ndim == 1:
        bc = structured_array_to_basic(coordinates_b).reshape((-1, 3))
    else:
        assert coordinates_b.ndim == 2
        bc = structured_array_to_basic(coordinates_b).reshape(
            (coordinates_b.shape[0], -1, 3))

    return coordinate_array_unaligned_broadcast_rmsd(ac, bc)


def coordinate_array_unaligned_broadcast_rmsd(coordinates_a, coordinates_b, out=None):
    """Calculate rmsd between entries in the given coordinate arrays without performing superposition.

    coordinates_a - array(shape=([a], n, 3), dtype=float)
    coordinates_b - array(shape=([b], n, 3), dtype=float)
    out - array(shape=(a, b), dtype=float)

    coordinates_a and coordinates_b must contain dimension n, the number of entries per coordinate set.

    returns broadcast_rmsd - array(shape=(a,b), dtype=float)
        Flattened to 1 dimension if a or b is dimension 2.
    """

    flatten = False

    if coordinates_a.ndim == 2:
        flatten = True
        coordinates_a = numpy.expand_dims(coordinates_a, 0)

    if coordinates_b.ndim == 2:
        flatten = True

        coordinates_b = numpy.expand_dims(coordinates_b, 0)

    if out is None:
        out = numpy.empty((coordinates_a.shape[0], coordinates_b.shape[0]))
    else:
        out = out.reshape((coordinates_a.shape[0], coordinates_b.shape[0]))

    out = numpy.sqrt(
        numpy.mean(
            numpy.square(
                numpy.linalg.norm(
                    numpy.expand_dims(coordinates_a, 1) -
                    numpy.expand_dims(coordinates_b, 0),
                    axis=-1)
            ),
            axis=-1)
    )

    if flatten:
        return out.ravel()
    else:
        return out
