// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>

//basic
#include <basic/Tracer.hh>

//numeric
#include <numeric/random/random.hh>

namespace devel {
namespace loop_creation {

static THREAD_LOCAL basic::Tracer TR( "devel.loop_creation.KICLoopCloser" );

//****CREATOR METHODS****//
std::string
KICLoopCloserCreator::keyname() const
{
	return KICLoopCloserCreator::mover_name();
}

protocols::moves::MoverOP
KICLoopCloserCreator::create_mover() const {
	return new KICLoopCloser;
}

std::string
KICLoopCloserCreator::mover_name()
{
	return "KICLoopCloser";
}

//****END CREATOR METHODS****//
/// @brief default constructor
KICLoopCloser::KICLoopCloser():
//	max_closure_attempts_(10),
//	max_KIC_moves_per_closure_attempt_(10000),
//	max_rama_score_increase_( 2.0 ),
//	max_total_delta_helix_( 15 ),
//	max_total_delta_strand_( 15 ),
//	max_total_delta_loop_( 15 ),
//	tolerance_( 0.01 )
{}

/// @brief explicit constructor
//KICLoopCloser::KICLoopCloser(
//	core::Size max_closure_attempts,
//	bool prevent_nonloop_modifications,
//	core::Size max_KIC_moves_per_closure_attempt,
//	core::Real max_rama_score_increase,
//	core::Real max_total_delta_helix,
//	core::Real max_total_delta_strand,
//	core::Real max_total_delta_loop,
//	core::Real tolerance
//):
//	LoopCloser(prevent_nonloop_modifications),
//	max_closure_attempts_(max_closure_attempts),
//	max_KIC_moves_per_closure_attempt_(max_KIC_moves_per_closure_attempt),
//	max_rama_score_increase_(max_rama_score_increase),
//	max_total_delta_helix_(max_total_delta_helix),
//	max_total_delta_strand_(max_total_delta_strand),
//	max_total_delta_loop_(max_total_delta_loop),
//	tolerance_(tolerance)
//{
//	init();
//}
	
protocols::moves::MoverOP
KICLoopCloser::clone() const {
	return( protocols::moves::MoverOP( new KICLoopCloser( *this ) ) );
}
protocols::moves::MoverOP
KICLoopCloser::fresh_instance() const {
	return protocols::moves::MoverOP( new KICLoopCloser );
}
	
std::string
KICLoopCloser::get_name() const {
	return "KICLoopCloser";
}
	
	
void
KICLoopCloser::apply(
	core::pose::Pose & pose
){
	success_=false;
	
	protocols::loops::loop_closure::kinematic_closure::KinematicMover kic_mover;
	protocols::loops::loop_closure::kinematic_closure::KinematicPerturberOP kic_perturber =
		new protocols::loops::loop_closure::kinematic_closure::NullKinematicPerturber(&kic_mover);
	kic_mover.set_perturber(kic_perturber);
	
	kic_mover.set_idealize_loop_first(true);
//	kic_mover->set_sample_nonpivot_torsions(false);

	core::Size n_pivot = loop().start();
	core::Size c_pivot = loop().stop();
	core::Size middle_pivot = numeric::random::random_range(n_pivot+1, c_pivot-1);
	TR << "Using following pivot residues for KIC: " <<  n_pivot << ", " <<
		middle_pivot << ", " << c_pivot << std::endl;
		
	TR << "n_pivot-1: " << pose.residue(n_pivot-1).atom("CA").xyz() << std::endl;
	TR << "n_pivot: " << pose.residue(n_pivot).atom("CA").xyz() << std::endl;
	TR << "c_pivot: " << pose.residue(c_pivot).atom("CA").xyz() << std::endl;
	TR << "c_pivot+1: " << pose.residue(c_pivot+1).atom("CA").xyz() << std::endl;
	TR << "middle_pivot: " << pose.residue(middle_pivot).atom("CA").xyz() << std::endl;

	core::conformation::Residue const & res=pose.residue(loop().cut());
	core::conformation::Residue const & next_res=pose.residue(loop().cut()+1);
	core::Real dist_squared = res.atom( res.upper_connect_atom() ).xyz().distance_squared(next_res.atom( next_res.lower_connect_atom() ).xyz());
	TR << "pre_KIC_broken distance: " << dist_squared << std::endl;

	//loop through all pivots here and check which one perturbs bb the least
	kic_mover.set_pivots(n_pivot, middle_pivot, c_pivot);
	kic_mover.apply(pose);

	if(kic_mover.last_move_succeeded())
	{
		TR << "n_pivot-1: " << pose.residue(n_pivot-1).atom("CA").xyz() << std::endl;
		TR << "n_pivot: " << pose.residue(n_pivot).atom("CA").xyz() << std::endl;
		TR << "c_pivot: " << pose.residue(c_pivot).atom("CA").xyz() << std::endl;
		TR << "c_pivot+1: " << pose.residue(c_pivot+1).atom("CA").xyz() << std::endl;
		TR << "middle_pivot: " << pose.residue(middle_pivot).atom("CA").xyz() << std::endl;
	
		core::conformation::Residue const & post_res=pose.residue(loop().cut());
		core::conformation::Residue const & post_next_res=pose.residue(loop().cut()+1);
		core::Real post_dist_squared = post_res.atom( post_res.upper_connect_atom() ).xyz().distance_squared(post_next_res.atom( post_next_res.lower_connect_atom() ).xyz());
		TR << "post_KIC_broken distance: " << post_dist_squared << std::endl;
		success_=true;
	}
}

/// @brief parse tag for use in RosettaScripts
void
KICLoopCloser::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose
){
	
}
	
} //loop_creation
} //devel
