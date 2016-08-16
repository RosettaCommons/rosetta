// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief simple mover for stealing side chains from one pose and sticking them
/// on another pose.
/// @author James Thompson

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/SequenceMapping.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/comparative_modeling/StealSideChainsMover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

StealSideChainsMover::StealSideChainsMover(
	core::pose::Pose const & source,
	core::id::SequenceMapping map
) :
	source_( source ),
	map_( map )
{}

void StealSideChainsMover::apply( core::pose::Pose & pose ) {
	using core::Size;
	using std::string;

	for ( Size ii = 1; ii <= map_.size1(); ++ii ) {
		Size const source_ii( map_[ii] );
		if ( source_ii == 0 ) continue;
		if ( source_ii > source_.total_residue() ) continue;

		string const name3( pose.residue_type(ii).name3() );
		string const name3_src( source_.residue_type(map_[ii]).name3());
		if ( name3 != name3_src ) continue;

		core::conformation::ResidueOP new_residue(
			source_.residue( map_[ii] ).clone()
		);

		// seqpos, new_residue, orient_backbone
		pose.replace_residue ( ii, *new_residue, true );
		//pose.replace_residue ( ii, *new_residue, false );
	}
} // apply

std::string
StealSideChainsMover::get_name() const {
	return "StealSideChainsMover";
}


} // comparative_modeling
} // protocols
