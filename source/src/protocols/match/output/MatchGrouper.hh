// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/MatchProcessor.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_MatchGrouper_hh
#define INCLUDED_protocols_match_output_MatchGrouper_hh

// Unit headers
#include <protocols/match/output/MatchGrouper.fwd.hh>

// Package headers
#include <protocols/match/Hit.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>

namespace protocols {
namespace match {
namespace output {

class MatchGrouper : public utility::pointer::ReferenceCount {
public:
	typedef core::Real Real;

public:
	MatchGrouper();

	virtual
	~MatchGrouper();

	virtual
	Size
	assign_group_for_match(
		match const & m
	) = 0;

	virtual
	Size
	assign_group_for_match(
		match_dspos1 const & m
	) = 0;

	virtual
	void
	reset() = 0;

};

}
}
}

#endif
