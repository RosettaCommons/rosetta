# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @author Sergey Lyskov

from __future__ import print_function

import rosetta, pyrosetta

pyrosetta.init()

if hasattr(rosetta, 'cereal'): # check if this is a cereliazation build

    pose = pyrosetta.pose_from_sequence("ARNDCEQGHILKMFPSTWYV")

    # saving into oss stream
    oss = rosetta.std.ostringstream()
    arc = rosetta.cereal.BinaryOutputArchive(oss);
    pose.save(arc)


    # loading from oss stream
    pose_l = pyrosetta.Pose()
    iss = rosetta.std.istringstream( oss.bytes() )
    arc = rosetta.cereal.BinaryInputArchive(iss);
    pose_l.load(arc)

    print(pose)
    print(pose_l)
    assert( pose_l.sequence() == pose.sequence() )
