#!/usr/bin/env python
# -*- coding: utf-8 -*-

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   GDock,py
## @brief  Determenistic Global Docking for Rosetta
## @author Sergey Lyskov


import rosetta
import rosetta.protocols.docking

import time

rosetta.init()


jumps = rosetta.utility.vector1_int()

pose = rosetta.pose_from_pdb('/Users/sergey/rosie_uploads/docking/1brs.pdb')
pose.pdb_info().name('D')  # changing name so PyMOL does not complain

scorefxn = rosetta.get_fa_scorefxn()
scorefxn(pose)


pymol = rosetta.PyMOL_Mover()
pymol.apply(pose)

rosetta.protocols.docking.setup_foldtree(pose, 'A_D', jumps)
docking_jump = jumps[1]


a1 = rosetta.numeric.xyzVector_double(0,0,0)
a2 = rosetta.numeric.xyzVector_double(1,0,0)
a3 = rosetta.numeric.xyzVector_double(0,1,0)
an = rosetta.numeric.xyzVector_double(0,0,1)


b1 = rosetta.numeric.xyzVector_double(10 + 0,0,0)
b2 = rosetta.numeric.xyzVector_double(10 + 1,0,0)
b3 = rosetta.numeric.xyzVector_double(10 + 0,1,0)
bn = rosetta.numeric.xyzVector_double(10 + 0,0,1)


# Checking if
# ((a2−a1 )×(a3−a1 )) ⋅ an ]] > 0
# ((b2−b1 )×(b3−b1 )) ⋅ bn ]] > 0  is hold

print ( a2 - a1 ).cross_product( a3 - a1 ).dot( an )
print ( b2 - b1 ).cross_product( b3 - b1 ).dot( bn )


# T(x) = Rx + s   where:
# R = U Q-1
# s = - U Q-1 a1 + b1

# Q = [a2-a1, a3-a1, an]
# U = [b2-b1, b3-b1, bn]


Q = rosetta.numeric.xyzMatrix_double()
U = rosetta.numeric.xyzMatrix_double()
#Q.row_x(a2-a1);  Q.row_y(a3-a1);  Q.row_z(an);
#U.row_x(b2-b1);  U.row_y(b3-b1);  U.row_z(bn);

Q.col_x(a2-a1);  Q.col_y(a3-a1);  Q.col_z(an);
U.col_x(b2-b1);  U.col_y(b3-b1);  U.col_z(bn);

R = U * rosetta.numeric.inverse(Q)
s = b1 - R * a1

def T(x): return R*x + s

print T(a1)
print T(a2)
print T(a3)

'''
for i in range(100):
    jump = pose.jump(docking_jump)
    jump.set_translation( rosetta.numeric.xyzVector_double(i,i,i) )
    pose.set_jump(docking_jump, jump)

    print scorefxn(pose)

    pymol.apply(pose)

    time.sleep(1)
'''
