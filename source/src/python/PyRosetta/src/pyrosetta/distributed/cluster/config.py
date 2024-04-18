# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.


__author__ = "Jason C. Klima"

import sys

from typing import List


environment_cmd: str = "conda env export --prefix {0}".format(
    sys.prefix
)  # YML file string generation command.
source_domains: List[str] = [
    "conda.graylab.jhu.edu",
    "west.rosettacommons.org",
    "conda.rosettacommons.org",
]  # Conda channels and/or source domains (containing PyRosetta usernames/passwords) to be stripped from YML file strings.
