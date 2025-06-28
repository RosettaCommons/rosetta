# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
# Bind warnings to methods


__author__ = "Jason C. Klima"


import functools
import numpy
import pyrosetta
import warnings

from pyrosetta.bindings.utility import bind_method
from pyrosetta.rosetta.core.pose import (
    setPoseExtraScore as _setPoseExtraScore,
    getPoseExtraFloatScores,
    getPoseExtraStringScores,
)


@bind_method(pyrosetta.rosetta.core.pose)
@functools.wraps(_setPoseExtraScore)
def setPoseExtraScore(pose=None, name=None, value=None):
    _setPoseExtraScore(pose, name, value)
    if isinstance(value, (float, int, bool)):
        result = getPoseExtraFloatScores(pose).__getitem__(name)
        if isinstance(value, float):
            diff = abs(result - value)
            try:
                numpy.testing.assert_equal(
                    result,
                    value,
                    err_msg=f"Pose extra score '{name}' float value precision is changing upon setting. DIFFERENCE: {diff}",
                    verbose=True,
                )
            except AssertionError as ex:
                warnings.warn(
                    str(ex).replace("\n", " "),
                    UserWarning,
                    stacklevel=2,
                )
        elif isinstance(value, (int, bool)):
            warnings.warn(
                f"Pose extra score '{name}' value is being type casted to float value: {result}",
                UserWarning,
                stacklevel=2,
            )
    elif isinstance(value, (bytes, bytearray)):
        try:
            result = getPoseExtraStringScores(pose).__getitem__(name)
        except UnicodeDecodeError as ex:
            warnings.warn(
                (
                    f"Pose extra score '{name}' value will be type casted to a `str` object that "
                    "cannot be decoded by UTF-8! This will result in a `UnicodeDecodeError` getting"
                    "raised when the value is accessed from PyRosetta. Consider serializing the value "
                    "into a `str` object and re-setting the serialized value to prevent the error."
                ),
                UserWarning,
                stacklevel=2,
            )
