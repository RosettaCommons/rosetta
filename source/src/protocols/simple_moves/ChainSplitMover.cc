// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/ChainSplitMover
/// @brief ChainSplitMover splits a pose into two chains given a cutpoint.
/// @author Robert Lindner <rlindner@mpimf-heidelberg.mpg.de>

// Unit Headers
#include <protocols/simple_moves/ChainSplitMover.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

namespace protocols {
namespace simple_moves {

ChainSplitMover::ChainSplitMover() :
	Mover( "ChainSplitMover" ),
	cutpoint_( 0 )
{}

ChainSplitMover::ChainSplitMover( core::Size cutpoint ) :
	Mover( "ChainSplitMover" ),
	cutpoint_( cutpoint )
{}

ChainSplitMover::ChainSplitMover( ChainSplitMover const & src ):
	Mover( "ChainSplitMover" ),
	cutpoint_( src.cutpoint_ )
{}

ChainSplitMover::~ChainSplitMover() = default;

void
ChainSplitMover::apply( core::pose::Pose & pose )
{
	// there's no point in setting a cutpoint before the first or after the last residue...
	assert( cutpoint_ > 0 && cutpoint_ < pose.total_residue() );

	pose.conformation().insert_chain_ending( cutpoint_ );
	core::kinematics::FoldTree jump_fold_tree( pose.fold_tree() );
	core::Size lower_half_middle = (cutpoint_ + 1) / 2;
	core::Size upper_half_middle = ( pose.total_residue()+ 1 - cutpoint_ )/2 + cutpoint_;
	jump_fold_tree.new_jump( lower_half_middle, upper_half_middle, cutpoint_ );
	pose.fold_tree( jump_fold_tree );
	core::pose::correctly_add_cutpoint_variants( pose );
}

protocols::moves::MoverOP
ChainSplitMover::clone() const
{
	return protocols::moves::MoverOP( new ChainSplitMover( *this ) );
}

protocols::moves::MoverOP
ChainSplitMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ChainSplitMover );
}

std::string
ChainSplitMover::get_name() const
{
	return "ChainSplitMover";
}

} // namespace simple_moves
} // namespace protocols

