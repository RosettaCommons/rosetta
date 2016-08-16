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
            "Rosetta predictions using the ab initio relax protocol"

        self.initialize_options_parser()
        self.parse_options(argv)
        self.setup()

#    def setup_input_data(self):
#        print "Checkout abrelax input dataset..."
#        try:
#            p = subprocess.Popen(
#                args=['svn', 'checkout', "--non-interactive", "--no-auth-cache",
#                      'https://svn.rosettacommons.org/source/trunk/mini.data/tests/scientific/cluster/features/sample_sources/abrelax/input/abrelax_lr5_input_structures_60_by_103_sequences_100908',
#                      'input/'],
#                stderr=subprocess.PIPE)
#            err = p.communicate()[1]
#            if err != '': print err
#        except OSError, e:
#            print >>sys.stderr, "Execution Failed:", e
#
#        self.mvars["decoy_input_silent_structures"] = os.path.join(
#            os.getcwd(), "input",
#            "abrelax_lr5_input_structures_60_by_103_sequences_100908.silent")


def main(argv):

    ss = SampleSource(argv)
    ss.submit()


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
