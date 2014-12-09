// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/MatchFilter.hh
/// @brief  Declaration for abstract class to filter matches.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_MatchFilter_hh
#define INCLUDED_protocols_match_output_MatchFilter_hh

// Unit headers
#include <protocols/match/output/MatchFilter.fwd.hh>

// Package headers
#include <protocols/match/Hit.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <iterator>
#include <string>

namespace protocols {
namespace match {
namespace output {

class MatchFilter : public utility::pointer::ReferenceCount {
public:

	MatchFilter( std::string filter_name );

	virtual
	~MatchFilter();

	/// @brief Returns true if the given match passes this filter
	virtual
	bool
	passes_filter(
		match const & m
	) const = 0;

	virtual
	bool
	passes_filter(
		match_dspos1 const & m
	) const = 0;

	std::string
	filter_name() const {
		return filter_name_; }

private:
	std::string filter_name_;

};

class StateAccumulatingMatchFilter : public MatchFilter {
public:

	StateAccumulatingMatchFilter( std::string filter_name );

	virtual
	~StateAccumulatingMatchFilter();

	/// @brief Returns true if the given match passes this filter
	virtual
	bool
	passes_filter(
		match const & m
	) const = 0;

	virtual
	bool
	passes_filter(
		match_dspos1 const & m
	) const = 0;

	/// @brief Note that a particular match has passed all the filters and will be output.
	virtual
	void
	note_match_accepted(
		match const & m
	) = 0;

	virtual
	bool
	note_match_accepted(
		match_dspos1 const & m
	) const = 0;

	/// @brief Erase all tracking data on which matches have already been output.
	virtual
	void
	reset() = 0;

};

}
}
}

#endif
