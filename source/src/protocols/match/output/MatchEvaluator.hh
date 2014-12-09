// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/MatchEvaluator.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_MatchEvaluator_hh
#define INCLUDED_protocols_match_output_MatchEvaluator_hh

// Unit headers
#include <protocols/match/output/MatchEvaluator.fwd.hh>

// Package headers
#include <protocols/match/Hit.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace match {
namespace output {

class MatchEvaluator : public utility::pointer::ReferenceCount
{
public:
	typedef core::Real Real;

public:
	virtual
	~MatchEvaluator();

	virtual
	Real
	score( match const & m ) const = 0;

	virtual
	Real
	score( match_dspos1 const & m ) const = 0;

};

}
}
}

#endif
