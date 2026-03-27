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

pyrosetta.init()

if hasattr(rosetta, 'cereal'): # check if this is a cereliazation build

    pose1 = pyrosetta.pose_from_sequence("ARNDCEQGHILKMFPSTWYV")
    pose2 = pyrosetta.pose_from_sequence("DDDDDD")

    # saving into oss stream
    oss = rosetta.std.ostringstream()
    arc = rosetta.cereal.BinaryOutputArchive(oss);
    pose1.save(arc)
    pose2.save(arc)


    # loading from oss stream
    pose1_l = pyrosetta.Pose()
    pose2_l = pyrosetta.Pose()
    iss = rosetta.std.istringstream( oss.bytes() )
    arc = rosetta.cereal.BinaryInputArchive(iss);
    pose1_l.load(arc)
    pose2_l.load(arc)

    print( 'Pose1  :\n{}-----------------\n'.format(pose1) )
    print( 'Pose1_l:\n{}-----------------\n'.format(pose1_l) )
    print( 'Pose2  :\n{}-----------------\n'.format(pose2) )
    print( 'Pose2_l:\n{}-----------------\n'.format(pose2_l) )
    assert( pose1_l.sequence() == pose1.sequence() )
    assert( pose2_l.sequence() == pose2.sequence() )
