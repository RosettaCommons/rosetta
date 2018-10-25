import numpy
import collections
import warnings
warnings.filterwarnings(
    "ignore",
    message="Numpy.*viewing or writing to an array returned by selecting multiple fields in a structured array.",
    category=FutureWarning,
    module='pyrosetta.utility.array'
)

def to_structured_array(*args, **kwargs):
    """Create structured array via (str->array) field mapping.

    Create structured array via concatenation of structured arrays or mapping of
    named fields.

    Args:
      *args:
        (array): Fields extracted from input structured array.
        (mapping(str, array)): Fields created from mapping entries.
        (tuple(str, array)): Single named field tuple

      *kwargs: Interpreted as mapping(str, array) argument.

      ndmin: Minimum dimension of result array.

    Returns:
      Structured array of extracted fields. Field values are broadcast.
    """

    ndmin = kwargs.pop("ndmin", 0)

    inputs = []
    args = list(args) + [kwargs]
    for a in args:
        if isinstance(a, numpy.ndarray):
            assert a.dtype.names is not None
            for n in a.dtype.names:
                inputs.append((n, numpy.array(a[n], copy=False, ndmin=ndmin)))
        elif isinstance(a, collections.Mapping):
            for n, v in a.items():
                assert isinstance(n, str)
                inputs.append((n, numpy.array(v, copy=False, ndmin=ndmin)))
        else:
            n, v = a
            assert isinstance(n, str)
            inputs.append((n, numpy.array(v, copy=False, ndmin=ndmin)))

    ndim = min(v.ndim for n, v in inputs)
    in_shape = [
        max(v.shape[d] for n, v in inputs) for d in range(ndim)  # noqa
    ]
    in_dtype = [
        (n, v.dtype, v.shape[ndim:]) for n, v in inputs  # noqa
    ]

    result = numpy.empty(in_shape, in_dtype)
    for n, v in inputs:
        result[n] = v

    return result


def extend_structured_array(source_array, extension_dtype):
    """Expand structured array with given fields.

    Returns structured array with fields of source_array and fields defined in
    extension_dtype. Fields present in source_array are initialzed with values of
    source_array. Fields in extension_dtype are zero initialized.
    """

    extension_dtype = numpy.dtype(extension_dtype)

    if not source_array.dtype.names:
        raise ValueError("source_array must be structured array")

    if not extension_dtype.names:
        raise ValueError("extension_dtype must be structured")

    intersecting_names = set(source_array.dtype.names).intersection(set(extension_dtype.names))
    if intersecting_names:
        raise ValueError("Extension contained existing names: %s" % intersecting_names)

    result = numpy.zeros_like(source_array, dtype=source_array.dtype.descr + extension_dtype.descr)

    rbuf = numpy.ndarray(
        shape=result.shape,
        dtype=source_array.dtype,
        buffer=result.data,
        strides=result.strides)
    rbuf[:] = source_array

    return result


def assign_by_field_names(to_array, from_array):
    """Copy structured array data by equivalent field names.

    Copy data from 'from_array' to 'to_array' by matching, potentially nested,
    field names. Partially replicates pre-1.13 numpy structured array
    assignment behavior.
    """

    common_names = set(from_array.dtype.names).intersection(
        set(to_array.dtype.names))

    structured_names = {
        n for n in common_names
        if from_array.dtype[n].names or to_array.dtype[n].names
    }

    for n in structured_names:
        assign_by_field_names(to_array[n], from_array[n])

    for n in common_names - structured_names:
        to_array[n] = from_array[n]


def structured_array_to_basic(in_array):
    """Convert structured array with fields of homogenous dtype to
    high-dimensional basic array.

    Convert structured array with fields of homogenous dtype to contiguous
    array of field dtype. Result array will be of shape ([in_shape],
    len(fields)) and of field dtype. Note, if field dtype is shaped, the field
    shape will be included as additional minor axis.

    Eg:
        in_array((10), [("a", "f", 3), ("b", "f", 3)]) ->
        out_array((10, 2, 3), "f")
    """

    if not in_array.dtype.fields:
        raise ValueError(
            "in_array must be structured", in_array.dtype)

    field_dtype = set(f[0] for f in in_array.dtype.fields.values())
    if len(field_dtype) > 1:
        raise ValueError(
            "in_array must be of homogenous field dtype", in_array.dtype)
    field_dtype = field_dtype.pop()

    if not field_dtype.subdtype:
        expected_dtype = field_dtype
        expected_shape = in_array.shape + (len(in_array.dtype.fields),)
    else:
        expected_dtype = field_dtype.subdtype[0]
        expected_shape = in_array.shape + (len(in_array.dtype.fields),) + field_dtype.subdtype[1]

    return numpy.ascontiguousarray(in_array).view(dtype=expected_dtype).reshape(expected_shape)


def basic_array_to_structured(in_array, field_names, field_dtype=None):
    """Convert basic array to structured array with fields of homogenous dtype."""

    if isinstance(field_names, str):
        field_names = (field_names,)

    if field_dtype is None:
        field_dtype = in_array.dtype
    else:
        field_dtype = numpy.dtype(field_dtype)

    result_dtype = numpy.dtype([(f, field_dtype) for f in field_names])
    if not field_dtype.subdtype:
        minor_shape = (len(field_names), )
    else:
        minor_shape = (len(field_names),) + field_dtype.subdtype[1]

    if in_array.shape[-len(minor_shape):] != minor_shape:
        raise ValueError(
            "in_array of invalid minor axis shape", in_array.shape[-len(minor_shape):])

    expected_shape = in_array.shape[:-len(minor_shape)]
    return in_array.reshape(expected_shape + (-1,)).view(result_dtype)


def atom_array_to_coordinates(atom_array):
    """Convert structured array of named atomic positions to coordinate buffer."""
    return structured_array_to_basic(atom_array).reshape(atom_array.shape[:-1] + (-1, 3))


def coordinate_array_to_atoms(coordinate_array, atom_names):
    """Convert coordinate buffer into structured array of named atomic positions."""
    return basic_array_to_structured(
        coordinate_array.reshape(coordinate_array.shape[:-2] + (-1, len(atom_names), 3)),
        field_names=atom_names, field_dtype=(coordinate_array.dtype, (3,))).squeeze(-1)


def rolling_window(array, window=(0,), asteps=None, wsteps=None, axes=None, toend=True):
    """Create a view of `array` which for every point gives the n-dimensional
    neighbourhood of size window. New dimensions are added at the end of
    `array` or after the corresponding original dimension.

    Parameters
    ----------
    array : array_like
        Array to which the rolling window is applied.
    window : int or tuple
        Either a single integer to create a window of only the last axis or a
        tuple to create it for the last len(window) axes. 0 can be used as a
        to ignore a dimension in the window.
    asteps : tuple
        Aligned at the last axis, new steps for the original array, ie. for
        creation of non-overlapping windows. (Equivalent to slicing result)
    wsteps : int or tuple (same size as window)
        steps for the added window dimensions. These can be 0 to repeat values
        along the axis.
    axes: int or tuple
        If given, must have the same size as window. In this case window is
        interpreted as the size in the dimension given by axes. IE. a window
        of (2, 1) is equivalent to window=2 and axis=-2.
    toend : bool
        If False, the new dimensions are right after the corresponding original
        dimension, instead of at the end of the array. Adding the new axes at the
        end makes it easier to get the neighborhood, however toend=False will give
        a more intuitive result if you view the whole array.

    Returns
    -------
    A view on `array` which is smaller to fit the windows and has windows added
    dimensions (0s not counting), ie. every point of `array` is an array of size
    window.

    Examples
    --------
    >>> a = numpy.arange(9).reshape(3,3)
    >>> rolling_window(a, (2,2))
    array([[[[0, 1],
             [3, 4]],

            [[1, 2],
             [4, 5]]],


           [[[3, 4],
             [6, 7]],

            [[4, 5],
             [7, 8]]]])

    Or to create non-overlapping windows, but only along the first dimension:
    >>> rolling_window(a, (2,0), asteps=(2,1))
    array([[[0, 3],
            [1, 4],
            [2, 5]]])

    Note that the 0 is discared, so that the output dimension is 3:
    >>> rolling_window(a, (2,0), asteps=(2,1)).shape
    (1, 3, 2)

    This is useful for example to calculate the maximum in all (overlapping)
    2x2 submatrixes:
    >>> rolling_window(a, (2,2)).max((2,3))
    array([[4, 5],
           [7, 8]])

    Or delay embedding (3D embedding with delay 2):
    >>> x = numpy.arange(10)
    >>> rolling_window(x, 3, wsteps=2)
    array([[0, 2, 4],
           [1, 3, 5],
           [2, 4, 6],
           [3, 5, 7],
           [4, 6, 8],
           [5, 7, 9]])
    """
    array = numpy.asarray(array)
    orig_shape = numpy.asarray(array.shape)
    window = numpy.atleast_1d(window).astype(int)

    if axes is not None:
        axes = numpy.atleast_1d(axes)
        w = numpy.zeros(array.ndim, dtype=int)
        for axis, size in zip(axes, window):
            w[axis] = size
        window = w

    # Check if window is legal:
    if window.ndim > 1:
        raise ValueError("`window` must be one-dimensional.")
    if numpy.any(window < 0):
        raise ValueError("All elements of `window` must be larger then 1.")
    if len(array.shape) < len(window):
        raise ValueError("`window` length must be less or equal `array` dimension.")

    _asteps = numpy.ones_like(orig_shape)
    if asteps is not None:
        asteps = numpy.atleast_1d(asteps)
        if asteps.ndim != 1:
            raise ValueError("`asteps` must be either a scalar or one dimensional.")
        if len(asteps) > array.ndim:
            raise ValueError("`asteps` cannot be longer then the `array` dimension.")
        # does not enforce alignment, so that steps can be same as window too.
        _asteps[-len(asteps):] = asteps

        if numpy.any(asteps < 1):
            raise ValueError("All elements of `asteps` must be larger then 1.")
    asteps = _asteps

    _wsteps = numpy.ones_like(window)
    if wsteps is not None:
        wsteps = numpy.atleast_1d(wsteps)
        if wsteps.shape != window.shape:
            raise ValueError("`wsteps` must have the same shape as `window`.")
        if numpy.any(wsteps < 0):
            raise ValueError("All elements of `wsteps` must be larger then 0.")

        _wsteps[:] = wsteps
        _wsteps[window == 0] = 1  # make sure that steps are 1 for non-existing dims.
    wsteps = _wsteps

    # Check that the window would not be larger then the original:
    if numpy.any(orig_shape[-len(window):] < window * wsteps):
        raise ValueError("`window` * `wsteps` larger then `array` in at least one dimension.")

    new_shape = orig_shape  # just renaming...

    # For calculating the new shape 0s must act like 1s:
    _window = window.copy()
    _window[_window == 0] = 1

    new_shape[-len(window):] += wsteps - _window * wsteps
    new_shape = (new_shape + asteps - 1) // asteps
    # make sure the new_shape is at least 1 in any "old" dimension (ie. steps
    # is (too) large, but we do not care.
    new_shape[new_shape < 1] = 1
    shape = new_shape

    strides = numpy.asarray(array.strides)
    strides *= asteps
    new_strides = array.strides[-len(window):] * wsteps

    # The full new shape and strides:
    if toend:
        new_shape = numpy.concatenate((shape, window))
        new_strides = numpy.concatenate((strides, new_strides))
    else:
        _ = numpy.zeros_like(shape)
        _[-len(window):] = window
        _window = _.copy()
        _[-len(window):] = new_strides
        _new_strides = _

        new_shape = numpy.zeros(len(shape) * 2, dtype=int)
        new_strides = numpy.zeros(len(shape) * 2, dtype=int)

        new_shape[::2] = shape
        new_strides[::2] = strides
        new_shape[1::2] = _window
        new_strides[1::2] = _new_strides

    new_strides = new_strides[new_shape != 0]
    new_shape = new_shape[new_shape != 0]

    return numpy.lib.stride_tricks.as_strided(array, shape=new_shape, strides=new_strides)
