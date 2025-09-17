# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.
# PyRosetta exceptions

__author__ = "Jason C. Klima"


class PyRosettaInitializationBaseException(RuntimeError):
    def __init__(self, message):
        super().__init__(message)
        self.message = message

    def __str__(self):
        return self.message

    def __repr__(self):
        return "{0}({1!r})".format(self.__class__.__name__, self.message)


class PyRosettaIsInitializedError(PyRosettaInitializationBaseException):
    def __init__(self, message):
        super().__init__(
            "{0}\nPlease ensure that PyRosetta is not already initialized (e.g., ".format(message)
            + "if using a Jupyter notebook, please restart the kernel) and try again."
        )


class PyRosettaIsNotInitializedError(PyRosettaInitializationBaseException):
    def __init__(self, message):
        super().__init__(
            "{0}\nPlease ensure that PyRosetta is initialized and try again.".format(message)
        )
