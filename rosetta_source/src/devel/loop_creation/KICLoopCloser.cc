// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file KICLoopCloser.cc
///
/// @brief Simple LoopCloser implementation that uses KIC to close the given loop. This closer
/// assumes that the loop is not currently closed and thus first idealizes from loop_start to loop_stop+1
/// @author Tim Jacobs

//unit
#include <devel/loop_creation/KICLoopCloser.hh>
#include <devel/loop_creation/KICLoopCloserCreator.hh>

//core
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

//TEMP
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/xyzVector.hh>

//protocols
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>

//basic
#include <basic/Tracer.hh>

//numeric
#include <numeric/random/random.hh>

namespace devel {
namespace loop_creation {

static basic::Tracer TR( "devel.loop_creation.KICLoopCloser" );

KICLoopCloserCreator::KICLoopCloserCreator() {}
KICLoopCloserCreator::~KICLoopCloserCreator() {}
LoopCloserOP KICLoopCloserCreator::create_loop_closer() const {
	return new KICLoopCloser;
}

std::string KICLoopCloserCreator::closer_name() const {
	return "KICLoopCloser";
}
	
KICLoopCloser::KICLoopCloser():
prevent_nonloop_changes_(true)
{}

bool
KICLoopCloser::close_loop(
	core::pose::Pose & pose,
	protocols::loops::Loop loop
){

	//prepare special foldtree
	core::kinematics::FoldTree saved_ft = pose.fold_tree();
	if(prevent_nonloop_changes_)
	{
		core::kinematics::FoldTree new_ft;
		new_ft.add_edge(1, loop.start()-1, core::kinematics::Edge::PEPTIDE);
		new_ft.add_edge(1, loop.start(), 1);
//		new_ft.add_edge(loop.start(), loop.stop()+1, core::kinematics::Edge::PEPTIDE);
//		new_ft.add_edge(1, loop.stop()+2, 2);
//		new_ft.add_edge(loop.stop()+2, pose.total_residue(), core::kinematics::Edge::PEPTIDE);
		new_ft.add_edge(loop.start(), pose.total_residue(), core::kinematics::Edge::PEPTIDE);
		if(!new_ft.check_fold_tree())
		{
			utility_exit_with_message("Bad fold tree created by KICLoopCloser. File a bug!");
		}
		
		TR.Debug << "Foldtree to prevent downstream propogation: " <<  new_ft << std::endl;
		pose.fold_tree(new_ft);
	}

	protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP kic_mover =
		new protocols::loops::loop_closure::kinematic_closure::KinematicMover();
	kic_mover->set_idealize_loop_first(false);
//	kic_mover->set_sample_nonpivot_torsions(false);

	core::Size n_pivot = loop.start();
	core::Size c_pivot = loop.stop()+1;
	core::Size middle_pivot = numeric::random::random_range(loop.start()+1, loop.stop()-1);
	TR << "Using following pivot residues for KIC: " <<  n_pivot << ", " <<
		middle_pivot << ", " << c_pivot << std::endl;
		
	TR << "n_pivot-1: " << pose.residue(n_pivot-1).atom("CA").xyz() << std::endl;
	TR << "n_pivot: " << pose.residue(n_pivot).atom("CA").xyz() << std::endl;
	TR << "c_pivot: " << pose.residue(c_pivot).atom("CA").xyz() << std::endl;
	TR << "c_pivot+1: " << pose.residue(c_pivot+1).atom("CA").xyz() << std::endl;
	TR << "middle_pivot: " << pose.residue(middle_pivot).atom("CA").xyz() << std::endl;

	core::conformation::Residue const & res=pose.residue(loop.cut());
	core::conformation::Residue const & next_res=pose.residue(loop.cut()+1);
	core::Real dist_squared = res.atom( res.upper_connect_atom() ).xyz().distance_squared(next_res.atom( next_res.lower_connect_atom() ).xyz());
	TR << "pre_KIC_broken distance: " << dist_squared << std::endl;

	//loop through all pivots here and check which one perturbs bb the least
	kic_mover->set_pivots(n_pivot, middle_pivot, c_pivot);
	kic_mover->apply(pose);

	//restore the saved fold tree
	pose.fold_tree(saved_ft);
	if(kic_mover->last_move_succeeded())
	{
		TR << "n_pivot-1: " << pose.residue(n_pivot-1).atom("CA").xyz() << std::endl;
		TR << "n_pivot: " << pose.residue(n_pivot).atom("CA").xyz() << std::endl;
		TR << "c_pivot: " << pose.residue(c_pivot).atom("CA").xyz() << std::endl;
		TR << "c_pivot+1: " << pose.residue(c_pivot+1).atom("CA").xyz() << std::endl;
		TR << "middle_pivot: " << pose.residue(middle_pivot).atom("CA").xyz() << std::endl;
	
		core::conformation::Residue const & post_res=pose.residue(loop.cut());
		core::conformation::Residue const & post_next_res=pose.residue(loop.cut()+1);
		core::Real post_dist_squared = post_res.atom( post_res.upper_connect_atom() ).xyz().distance_squared(post_next_res.atom( post_next_res.lower_connect_atom() ).xyz());
		TR << "post_KIC_broken distance: " << post_dist_squared << std::endl;
		return true;
	}
	return false;
}

///@brief parse tag for use in RosettaScripts
void
KICLoopCloser::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose
){
	
}
	
} //loop_creation
} //devel
