// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/UpstreamBuilder.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_PoseInserter_hh
#define INCLUDED_protocols_match_output_PoseInserter_hh

// Unit headers
#include <protocols/match/output/PoseInserter.fwd.hh>

// Package headers
#include <protocols/match/upstream/UpstreamBuilder.hh>

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {

class PoseInserter : public upstream::UpstreamResidueProcessor
{
public:
	typedef core::pose::Pose Pose;
	typedef core::Size       Size;
public:
	PoseInserter( Pose & pose_to_modify );
	PoseInserter( Pose & pose_to_modify, Size resid_to_replace );
	virtual ~PoseInserter();

	/// @brief Take a conformation::Residue from the upstream builder and
	/// call Pose::replace_residue at a particular position.
	virtual
	void
	process_hit(
		Hit const & hit,
		core::conformation::Residue const & upstream_conformation
	);

	void
	set_replacement_resid( Size seqpos );

private:
	Pose & pose_;
	Size   resid_to_replace_;
};

}
}
}

#endif
