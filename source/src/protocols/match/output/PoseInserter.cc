// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/PoseInserter.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/PoseInserter.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/exit.hh>

#include <utility/vector1.hh>


// C++ headers

namespace protocols {
namespace match {
namespace output {

PoseInserter::PoseInserter( Pose & pose_to_modify ) :
	pose_( pose_to_modify ),
	resid_to_replace_( 0 )
{}

PoseInserter::PoseInserter( Pose & pose_to_modify, Size resid_to_replace ) :
	pose_( pose_to_modify ),
	resid_to_replace_( resid_to_replace )
{}

PoseInserter::~PoseInserter() = default;

/// @brief Take a conformation::Residue from the upstream builder and
/// call Pose::replace_residue at a particular position.
void
PoseInserter::process_hit(
	Hit const &,
	core::conformation::Residue const & upstream_conformation
)
{
	runtime_assert( resid_to_replace_ != 0 );

	pose_.replace_residue( resid_to_replace_, upstream_conformation, false );

}

void
PoseInserter::set_replacement_resid( Size seqpos )
{
	resid_to_replace_ = seqpos;
}

}
}
}

