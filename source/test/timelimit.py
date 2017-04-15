#!/usr/bin/env python
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   timelimit.py
## @brief  shell wraper to limit execution time
## @author Sergey Lyskov
from __future__ import print_function

import sys
if not hasattr(sys, "version_info") or sys.version_info < (2,4):
    raise ValueError("Script requires Python 2.4 or higher!")

import subprocess, signal, time, sys, os

def main(argv):
    #print argv
    if len(argv) < 3:
        print('Usage: "timelimit <number of minutes to wait> <command line> <args...>"')
        return 1

    timeout = int(argv[1])

    commline = ' '.join(argv[2:])
    proc = subprocess.Popen(["bash", "-c", commline])
    start = time.time()
    while time.time() - start <= timeout*60:
        retcode = proc.poll()
        if retcode is not None: break
        time.sleep(5)
    if retcode is None:
        print("*** '%s' exceeded the timeout and will be killed!" % commline)
        os.kill(proc.pid, signal.SIGTERM)
        return 1
    return retcode


if __name__ == "__main__":
    sys.exit(main(sys.argv))
