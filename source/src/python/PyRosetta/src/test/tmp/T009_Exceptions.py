# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

import pyrosetta
import pyrosetta.rosetta as rosetta

import sys

pyrosetta.init()

xml = "<ROSETTASCRIPTS></ROSETTASCRIPTS>"
tag = rosetta.utility.tag.Tag.create(xml)

options = rosetta.basic.options.process()
parser = rosetta.protocols.rosetta_scripts.RosettaScriptsParser()
pose = pyrosetta.Pose()

try:
    parser.generate_mover_for_protocol(pose, False, tag, options)
except:
    t = sys.exc_info()[0]
    e = sys.exc_info()[1]
    print('Exceptions type: {}, message: {}'.format(t, e) )
    assert 'parser::protocol file must specify PROTOCOLS section' in str(e)

try:
    parser.generate_mover_for_protocol(pose, False, tag, options)
except RuntimeError as e:
    assert 'parser::protocol file must specify PROTOCOLS section' in str(e)
