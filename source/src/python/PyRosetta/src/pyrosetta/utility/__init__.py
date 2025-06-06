# :noTabs=true:
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington CoMotion, email: license@uw.edu.

__author__ = "Jason C. Klima"


def get_package_version(distribution_name):
    from functools import singledispatch
    from importlib.metadata import PackageNotFoundError, version

    @singledispatch
    def to_tuple(v):
        raise ValueError(
            f"Version type {type(v)} not supported for '{distribution_name}' version: {v}"
        )

    @to_tuple.register(str)
    def _from_str(v):
        try:
            t = tuple(map(int, v.split(".")[:3]))
        except ValueError:
            t = (v,)
        finally:
            return t

    @to_tuple.register(int)
    def _from_int(v):
        return (v,)

    @to_tuple.register(type(None))
    def _from_none(none):
        return None

    try:
        _version = version(distribution_name)
    except PackageNotFoundError:
        _version = None

    return to_tuple(_version)
