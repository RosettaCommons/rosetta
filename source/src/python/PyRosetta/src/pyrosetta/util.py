# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# Developer utility functions
# @author atom-moyer

def bind_method(class_):
    """
    bind_method is to be used as a decorator for binding methods
    onto classes which are defined in C++. This functionality is
    necessary because C++ objects are not explicitly defined in the python code.
    Therefore, dynamically defining is the only option.  Binding methods
    should be limited to python specific challenges.

    @atom-moyer
    """
    def _bind_method(method):
        setattr(class_, method.__name__, method)
        return method
    return _bind_method


def vector1_indices(slice_, size):
    """
    takes an instance of the built-in slice object and a length of a slice-able
    object and returns index 1 compatible indices. Think of this as
    `slice.vector1_indices(size)` which is analogous to `slice.indices(size)`.

    @atom-moyer
    """
    return slice(slice_.start if slice_.start else 1,
                 slice_.stop if slice_.stop else size + 1,
                 slice_.step).indices(size + 1)
