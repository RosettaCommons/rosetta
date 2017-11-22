# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

print('-------- Test/Demo for capturing Tracers output in PyRosetta --------')

import rosetta, pyrosetta


rosetta.basic.Tracer.super_mute(False)

T = rosetta.basic.PyTracer()
rosetta.basic.Tracer.set_ios_hook(T, rosetta.basic.Tracer.get_all_channels_string(), False)

pyrosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')


pose = rosetta.core.import_pose.pose_from_file("../test/data/test_in.pdb")

print( '\nCaptured IO:' )
print( T.buf() )


# More fancy example, using a output callback:

class MyPyTracer(rosetta.basic.PyTracer):
    def __init__(self):
        rosetta.basic.PyTracer.__init__(self)
        self.set_ios_hook = rosetta.basic.Tracer.set_ios_hook

    def __del__(self):
        self.set_ios_hook(None, '')

    def output_callback(self, s):
        print('MyPyTracer.output_callback with argument:')
        print(s)

M = MyPyTracer()
rosetta.basic.Tracer.set_ios_hook(M, rosetta.basic.Tracer.get_all_channels_string())

pyrosetta.init()
