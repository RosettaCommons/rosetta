// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/components/NullPoseFolder.cc
/// @brief Pose folder that does nothing
/// @author Tom Linsky (tlinsky@uw.edu)

#include <protocols/denovo_design/components/NullPoseFolder.hh>

// Protocol headers

// Core headers

// Basic/Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.NullPoseFolder" );

namespace protocols {
namespace denovo_design {
namespace components {

NullPoseFolder::NullPoseFolder():
	PoseFolder( NullPoseFolder::class_name() )
{
}

NullPoseFolder::~NullPoseFolder()
{}

NullPoseFolder::PoseFolderOP
NullPoseFolder::clone() const
{
	return NullPoseFolder::PoseFolderOP( new NullPoseFolder( *this ) );
}

std::string
NullPoseFolder::class_name()
{
	return "NullPoseFolder";
}

void
NullPoseFolder::parse_tag( utility::tag::TagCOP , basic::datacache::DataMap & )
{
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
NullPoseFolder::apply(
	core::pose::Pose & ,
	core::select::residue_selector::ResidueSubset const & ,
	protocols::loops::Loops const & ) const
{
}

} //protocols
} //denovo_design
} //components
