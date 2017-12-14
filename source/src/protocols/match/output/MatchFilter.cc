// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/MatchFilter.cc
/// @brief  Implementation for abstract filter class
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/MatchFilter.hh>

// Package headers
#include <protocols/match/Hit.fwd.hh>

// Utility headers
#include <utility>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace match {
namespace output {

MatchFilter::MatchFilter( std::string const & filter_name )
: filter_name_(filter_name)
{}

MatchFilter::~MatchFilter() = default;

StateAccumulatingMatchFilter::StateAccumulatingMatchFilter( std::string filter_name )
: MatchFilter( filter_name )
{}

StateAccumulatingMatchFilter::~StateAccumulatingMatchFilter() = default;

}
}
}
