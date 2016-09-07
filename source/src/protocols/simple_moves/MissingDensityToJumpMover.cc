// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/MissingDensityToJumpMover.cc
/// @brief  Implementation of mover that inserts a jump where there is gap in the pdb. This gap corresponds to missing density.
/// @author TJ Brunette (tjbrunette@gmail.com), May 2011

// Unit Headers
#include <protocols/simple_moves/MissingDensityToJumpMover.hh>

// P
// tracer
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Residue.hh>

#include <core/kinematics/FoldTree.hh>

#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace simple_moves {
static THREAD_LOCAL basic::Tracer TR( "protocols.mover.MissingDensityToJumpMover" );

// Default constructor
MissingDensityToJumpMover::MissingDensityToJumpMover(): protocols::moves::Mover( MissingDensityToJumpMover::get_name() )
{}

// Copy constructor
MissingDensityToJumpMover::MissingDensityToJumpMover(MissingDensityToJumpMover const & ) = default;

MissingDensityToJumpMover::~MissingDensityToJumpMover() = default;

void
MissingDensityToJumpMover::apply( core::pose::Pose & pose ) {
	using namespace core;
	using namespace core::conformation;
	using namespace core::chemical;
	Size const nres( pose.total_residue() );
	for ( Size i=1; i< nres; ++i ) {//don't have to go to last residue thus < rather than <=
		if ( pose.residue_type(i).is_polymer() && !pose.residue_type(i).is_lower_terminus() && !pose.fold_tree().is_cutpoint(i) ) {
			Residue const &current_rsd(pose.residue(i));
			Residue const &next_rsd(pose.residue(i+1));
			core::Real bondlength = ( current_rsd.atom( current_rsd.upper_connect_atom() ).xyz() - next_rsd.atom( next_rsd.lower_connect_atom() ).xyz() ).length();
			if ( bondlength > 2.5 ) {
				TR << "[ WARNING ] missing density found at residue " << i << std::endl;
				core::kinematics::FoldTree update_tree(pose.fold_tree());
				update_tree.new_jump(i,i+1,i);
				pose.fold_tree(update_tree);
				core::pose::add_variant_type_to_pose_residue( pose, chemical::UPPER_TERMINUS_VARIANT, i);
				core::pose::add_variant_type_to_pose_residue( pose, chemical::LOWER_TERMINUS_VARIANT, i+1 );
			}
		}
	}
}

std::string MissingDensityToJumpMover::get_name() const {
	return "MissingDensityToJumpMover";
}

protocols::moves::MoverOP
MissingDensityToJumpMover::clone() const
{
	return protocols::moves::MoverOP( new MissingDensityToJumpMover(*this) );
}

protocols::moves::MoverOP
MissingDensityToJumpMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new MissingDensityToJumpMover() );
}

} // moves
} // protocols
