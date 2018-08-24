# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
"""Utility functions for method and property bindings over wrapper classes."""



def bind_method(class_):
    """
    bind_method is to be used as a decorator for binding methods
    onto classes which are defined in C++. This functionality is
    necessary because C++ objects are not explicitly defined in the python code.
    Therefore, dynamically defining is the only option.  Binding methods
    should be limited to python specific challenges.
    """
    def _bind_method(method):
        setattr(class_, method.__name__, method)
        return method
    return _bind_method


def bind_property(class_):
    """
    bind_property is to be used as a decorator for binding properties
    onto classes which are defined in C++. This functionality is
    necessary because C++ objects are not explicitly defined in the python code.
    Therefore, dynamically defining is the only option.  Binding properties
    should be limited to python specific challenges.
    """
    def _bind_property(method):
        setattr(class_, method.__name__, property(method))
        return method
    return _bind_property


def bind_classmethod(class_):
    """
    bind_classmethod is to be used as a decorator for binding class methods
    onto classes which are defined in C++. This functionality is
    necessary because C++ objects are not explicitly defined in the python code.
    Therefore, monkey patching is the only option.  Binding class methods
    should be limited to python specific challenges.
    """
    def _bind_classmethod(method):
        setattr(class_, method.__name__, classmethod(method))
        return method
    return _bind_classmethod


def slice_1base_indicies(slice_, size):
    """Convert a slice to 1-based indicies.
    takes an instance of the built-in slice object and a length of a slice-able
    object and returns index 1 compatible indices. Think of this as
    `slice.vector1_indices(size)` which is analogous to `slice.indices(size)`.
    """

    if slice_.start == 0:
        raise ValueError(
            "1-based slicing does not support 0 start indices: %s" % slice_)

    return slice(slice_.start if slice_.start else 1, slice_.stop, slice_.step).indices(size + 1)
