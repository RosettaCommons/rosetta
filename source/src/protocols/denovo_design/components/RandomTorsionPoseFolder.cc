// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/RandomTorsionPoseFolder.cc
/// @brief Folds a pose using random phi/psi torsions
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/components/RandomTorsionPoseFolder.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// Basic/Utililty headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>

static basic::Tracer TR( "protocols.denovo_design.components.RandomTorsionPoseFolder" );

namespace protocols {
namespace denovo_design {
namespace components {

RandomTorsionPoseFolder::RandomTorsionPoseFolder():
	PoseFolder( RandomTorsionPoseFolder::class_name() )
{
}

RandomTorsionPoseFolder::~RandomTorsionPoseFolder() = default;

RandomTorsionPoseFolder::PoseFolderOP
RandomTorsionPoseFolder::clone() const
{
	return RandomTorsionPoseFolder::PoseFolderOP( new RandomTorsionPoseFolder( *this ) );
}

std::string
RandomTorsionPoseFolder::class_name()
{
	return "RandomTorsionPoseFolder";
}

void
RandomTorsionPoseFolder::parse_tag( utility::tag::TagCOP , basic::datacache::DataMap & )
{
}

/// @brief checks to ensure a residue is in the ResidueSubset
///        exits with error if not
void
check_residue(
	core::select::residue_selector::ResidueSubset const & subset,
	core::Size const resid )
{
	if ( ( resid < 1 ) || ( resid > subset.size() ) ) {
		std::stringstream msg;
		msg << RandomTorsionPoseFolder::class_name() << "::apply(): A residue in a given loop (" << resid
			<< ") does not exist in the movable residue subset (" << subset.size() << " residues)" << std::endl;
		utility_exit_with_message( msg.str() );
	}
}

/// @brief checks to ensure loops are contained within the ResidueSubset
///        exits with error if not
void
check_residue_subset(
	core::select::residue_selector::ResidueSubset const & subset,
	protocols::loops::Loops const & loops )
{
	for ( auto const & loop : loops ) {
		check_residue( subset, loop.start() );
		check_residue( subset, loop.stop() );
	}
}

/// @brief performs folding
/// @param pose    - The pose to be folded, with all residues added.  The pose should be prepared with
///                  any necessary cutpoints added before giving to the PoseFolder. Torsions in the pose
///                  should be adjusted, and no residues should be added or removed.
/// @param movable - Subset of residues for which new backbone conformations will be sampled. Residues
///                  specified as 'True' in movable must also be present in one or more Loops in order
///                  to be folded. Movable's size must match pose.size()
/// @param loops   - Loops to be folded.  Cutpoints specified here must be match the cutpoints found in
///                  the pose. Residues not within any loop should not be folded. Residues contained
///                  in a loop but not in the movable set should not be folded.
/// @throws EXCN_Fold if anything goes wrong in folding. Derived classes should throw this.
void
RandomTorsionPoseFolder::apply(
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & movable,
	protocols::loops::Loops const & loops ) const
{
	check_residue_subset( movable, loops );

	for ( auto const & loop : loops ) {
		for ( core::Size resid=loop.start(); resid<=loop.stop(); ++resid ) {
			if ( !movable[resid] ) continue;
			TR.Debug << "Setting random torsions for residue << " << resid << std::endl;
			pose.set_phi( resid, numeric::random::rg().uniform()*360.0 - 180.0 );
			pose.set_psi( resid, numeric::random::rg().uniform()*360.0 - 180.0 );
		}
	}

	// remodel doesn't remove cutpoints, so this shouldn't either
	/*
	for ( protocols::loops::Loops::const_iterator l=loops.begin(); l!=loops.end(); ++l ) {
	if ( !l->cut() ) continue;
	core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, l->cut() );
	core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, l->cut() + 1 );
	}
	*/
}

} //protocols
} //denovo_design
} //components
