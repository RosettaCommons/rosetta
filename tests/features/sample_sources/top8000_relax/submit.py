#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os, sys
import subprocess

sys.path.insert(0, "/".join(os.getcwd().split("/")[:-1]))

from submit import BaseSampleSource

class SampleSource(BaseSampleSource):


    def __init__(self, argv):
        BaseSampleSource.__init__(self)

        self.sample_source_description_default_value = \
            "Top 8000 chains in PDB from the Richardson Lab with Reduce placed Hydrogens, relaxed with fast relax"

        self.initialize_options_parser()
        self.parse_options(argv)
        self.setup()


def main(argv):

    ss = SampleSource(argv)
    ss.submit()


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
