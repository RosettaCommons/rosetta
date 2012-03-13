// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody2/Ab_util.cc
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)

#include <protocols/antibody2/Ab_util.hh>

// Rosetta Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/id/AtomID_Map.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/import_pose/import_pose.hh>

#include <iostream>
#include <fstream>



static basic::Tracer TR("antibody2.Ab_util");


using namespace core;
namespace protocols{
namespace antibody2{

    
    //JQX:
    // OK Description (Jason contributed)
    // assuming a loop 10,11,12,13,14,15, cut_point is 13/14
    //                  or loop = Loop (10, 15, 13)
    // 1. (Rosetta Default) 
    //           set_single_loop_fold_tree()
    //           jump from 8->17 (-2 amd +2)
    // 2. (Aroop's) 
    //            simple_one_loop_fold_tree() 
    //            jump from 9->16 (-1 and +1)
    // 3. (Aroop's) 
    //            simple_fold_tree() 
    //            jump from 10->15, but you need manuly input
    // 4. (Aroop's)
    //            setup_simple_fold_tree()
    //            jump from 10->15, but you need manuly input
    
    void simple_one_loop_fold_tree(
                                   pose::Pose & pose_in,
                                   loops::Loop const & loop	) {
        using namespace kinematics;
        
        TR <<  "H3M Setting up simple one loop fold tree" << std::endl;
        
        //setup fold tree for this loop
        FoldTree f;
        f.clear();
        Size nres = pose_in.total_residue();
        Size jumppoint1 = loop.start() - 1;
        Size jumppoint2 = loop.stop() + 1;
        
        if( jumppoint1 < 1 )   jumppoint1 = 1;
        if( jumppoint2 > nres) jumppoint2 = nres;
        
        f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
        f.add_edge( jumppoint1, loop.cut(), Edge::PEPTIDE );
        f.add_edge( loop.cut() + 1, jumppoint2, Edge::PEPTIDE );
        f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
        f.add_edge( jumppoint1, jumppoint2, 1 );
        f.reorder( 1 );
        
        pose_in.fold_tree( f );
        
        TR <<  "H3M Finished setting up simple one loop fold tree" << std::endl;
        
        return;
    } // simple_one_loop_fold_tree
    
    
    
    
    
    void simple_fold_tree(
                          pose::Pose & pose_in,
                          Size jumppoint1,
                          Size cutpoint,
                          Size jumppoint2 ) {
        using namespace kinematics;
        
        TR <<  "H3M Setting up simple fold tree" << std::endl;
        
        //setup fold tree for this loop
        FoldTree f;
        f.clear();
        Size nres = pose_in.total_residue();
        
        if( jumppoint1 < 1 )   jumppoint1 = 1;
        if( jumppoint2 > nres) jumppoint2 = nres;
        
        f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
        f.add_edge( jumppoint1, cutpoint, Edge::PEPTIDE );
        f.add_edge( cutpoint + 1, jumppoint2, Edge::PEPTIDE );
        f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
        f.add_edge( jumppoint1, jumppoint2, 1 );
        f.reorder( 1 );
        
        pose_in.fold_tree( f );
        
        TR <<  "H3M Finished setting up simple fold tree" << std::endl;
        
        return;
    } // simple_fold_tree
    
    
    
    
    
    
	void setup_simple_fold_tree(
                                          Size jumppoint1,
                                          Size cutpoint,
                                          Size jumppoint2,
                                          Size nres,
                                          pose::Pose & pose_in ) {
        
		using namespace kinematics;
        
		TR << "ABM Setting up simple fold tree" << std::endl;
        
		FoldTree f;
		f.clear();
        
		f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
		f.add_edge( jumppoint1, cutpoint, Edge::PEPTIDE );
		f.add_edge( cutpoint + 1, jumppoint2, Edge::PEPTIDE );
		f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
		f.add_edge( jumppoint1, jumppoint2, 1 );
		f.reorder( 1 );
        
		pose_in.fold_tree( f );
        
		TR << "ABM Done: Setting up simple fold tree" << std::endl;
        
	} // setup_simple_fold_tree
    
    
    


} // namespace antibody2
} // namespace protocols






