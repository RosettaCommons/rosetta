// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief      Optimizes the protein embedding in the membrane
/// @details Optimizes the protein embedding in the membrane given the smooth
///   high-res score function; transforms the protein into the membrane,
///   optimizes the membrane position (flexible), and uses the optimized
///   embedding to reposition the protein in the membrane
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_OptimizeProteinEmbeddingMover_cc
#define INCLUDED_protocols_membrane_OptimizeProteinEmbeddingMover_cc

// Unit Headers
#include <protocols/membrane/OptimizeProteinEmbeddingMover.hh>
#include <protocols/membrane/OptimizeProteinEmbeddingMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <protocols/membrane/OptimizeMembranePositionMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <protocols/membrane/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.membrane.OptimizeProteinEmbeddingMover" );

namespace protocols {
namespace membrane {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Defaults: scorefxn = smooth2012
OptimizeProteinEmbeddingMover::OptimizeProteinEmbeddingMover() : protocols::moves::Mover()
{
	register_options();
}

/// @brief Destructor
OptimizeProteinEmbeddingMover::~OptimizeProteinEmbeddingMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
OptimizeProteinEmbeddingMover::clone() const {
	return ( protocols::moves::MoverOP( new OptimizeProteinEmbeddingMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
OptimizeProteinEmbeddingMover::fresh_instance() const {
	return protocols::moves::MoverOP( new OptimizeProteinEmbeddingMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
OptimizeProteinEmbeddingMover::parse_my_tag(
	utility::tag::TagCOP /*tag*/,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// TODO: implement this

}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
OptimizeProteinEmbeddingMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new OptimizeProteinEmbeddingMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
OptimizeProteinEmbeddingMoverCreator::keyname() const {
	return OptimizeProteinEmbeddingMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
OptimizeProteinEmbeddingMoverCreator::mover_name() {
	return "OptimizeProteinEmbeddingMover";
}


/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (OptimizeProteinEmbeddingMover)
std::string
OptimizeProteinEmbeddingMover::get_name() const {
	return "OptimizeProteinEmbeddingMover";
}

/// @brief Flip the downstream partner in the membrane
void OptimizeProteinEmbeddingMover::apply( Pose & pose ) {

	using namespace numeric;
	using namespace core::conformation::membrane;
	using namespace protocols::membrane;
	using namespace protocols::membrane::geometry;

	TR << "Start OptimizeProteinEmbeddingMover, transforming into membrane" << std::endl;
	TR << "WARNING: the membrane information in the MEM residue is lost!!!" << std::endl;

	// if pose not a membrane protein, add membrane
	if ( ! pose.conformation().is_membrane() ) {

		TR << "WARNING: You are calling the OptimizeProteinEmbeddingMover on a non-membrane pose. Creating a membrane protein from the pose..." << std::endl;
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );
	}

	// remember foldtree
	core::kinematics::FoldTree orig_ft = pose.fold_tree();

	// if membrane at root, reorder foldtree to have pose TM COM at root
	if ( is_membrane_fixed( pose ) ) {
		TR << "Reordering foldtree:" << std::endl;
		core::Size anchor( create_membrane_foldtree_anchor_pose_tmcom( pose ) );
		core::kinematics::FoldTree ft = pose.fold_tree();
		ft.reorder( anchor );
		pose.fold_tree( ft );
	}

	// starting foldtree
	TR << "Starting foldtree: Is membrane fixed? " << is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

	// transforms into fixed membrane and then optimizes MEM (=dangling)
	OptimizeMembranePositionMoverOP opt( new OptimizeMembranePositionMover() );
	opt->apply( pose );

	// get the optimized embedding from flexible MEM
	core::Vector center( pose.conformation().membrane_info()->membrane_center() );
	core::Vector normal( pose.conformation().membrane_info()->membrane_normal() );
	EmbeddingDefOP current_emb( new EmbeddingDef( center, normal ) );

	// use current embedding and transform into default membrane (fixed)
	TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover( current_emb ) );
	transform->use_default_membrane( true );
	transform->apply( pose );

	// reset the membrane residue because it was changed during optimization
	SetMembranePositionMoverOP set( new SetMembranePositionMover() );
	set->apply( pose );

	// print tilt angle and distance from membrane center
	pose_tilt_angle_and_center_distance( pose );

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree after reset: Is membrane fixed? " << is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

}// apply

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
void OptimizeProteinEmbeddingMover::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::mp::transform::optimize_embedding );

}

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_OptimizeProteinEmbeddingMover_cc
