# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @author Sergey Lyskov

from rosetta import *
rosetta.init(extra_options = "-constant_seed")  # WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION RUNS!!!!!
import os; os.chdir('.test.output')

print 'Fold tree --------------------------------------------------'

pose = pose_from_pdb("../test/data/test_in.pdb")


#NEEDED: Fold tree, getting, printing, setting jumps and cuts
#fold tree functions
existing_ft = pose.fold_tree()
print existing_ft

# I'm guessing on the below!!
ft = FoldTree()
ft.simple_tree(116)
ft.new_jump(30,50,40) # jump_begin, jump_end, cutpoint
print ft # FOLD_TREE  EDGE 1 30 -1  EDGE 30 40 -1  EDGE 30 50 1  EDGE 50 116 -1  EDGE 50 41 -1
ft.check_fold_tree()
print ft.nres(), ft.size(), ft.root() # 116, 5, 1

ft.clear()
ft.add_edge(1,30,-1)
ft.add_edge(30,40,-1)
ft.add_edge(30,50,1)
ft.add_edge(50,41,-1)
ft.add_edge(50,116,-1)
ft.check_fold_tree()
print ft

ft.clear()
# This function no longer present in C++
# ft.simple_fold_tree()
# ft.add_jump(30,50,40)
print ft
