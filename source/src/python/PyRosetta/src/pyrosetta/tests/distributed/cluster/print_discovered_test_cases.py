# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

import argparse
import os
import unittest

def flatten(suite):
    for test in suite:
        yield from (flatten(test) if isinstance(test, unittest.TestSuite) else [test])

def main(root):
    suite = unittest.defaultTestLoader.discover(root)
    test_cases = (test_case.id() for test_case in flatten(suite))
    print(*test_cases, sep="\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=str, default=os.curdir, help="Root directory path (Default: '.')")
    args = parser.parse_args()
    main(args.root)
