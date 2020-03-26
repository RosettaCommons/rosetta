// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/MatchFilter.hh
/// @brief  Declaration for abstract class to filter matches.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_LimitHitsPerRotamerFilter_hh
#define INCLUDED_protocols_match_output_LimitHitsPerRotamerFilter_hh

// Unit headers
#include <protocols/match/output/LimitHitsPerRotamerFilter.fwd.hh>

// Package headers
#include <protocols/match/output/MatchFilter.hh>
#include <protocols/match/Hit.fwd.hh>

// Utility headers

// C++ headers
#include <map>

#include <utility/OrderedTuple.fwd.hh>

#ifdef PYROSETTA
#include <utility/OrderedTuple.hh>
#endif

namespace protocols {
namespace match {
namespace output {

class LimitHitsPerRotamerFilter : public StateAccumulatingMatchFilter {
public:
	typedef core::Size Size;
	typedef std::map< utility::OrderedTuple< utility::vector1< core::Size > >, core::Size > RotamerComboCountMap;

public:
	LimitHitsPerRotamerFilter();
	LimitHitsPerRotamerFilter( core::Size n_geometric_constraints );

	void
	set_n_geometric_constraints( core::Size n_csts );

	void
	set_limit_for_rotamer_combo( core::Size limit );

	~LimitHitsPerRotamerFilter() override;

	/// @brief Returns true if the given match passes this filter
	bool
	passes_filter(
		match const & m
	) const override;

	/// @brief Note that a particular match has passed all the filters and will be output.
	void
	note_match_accepted(
		match const & m
	) override;

	/// @brief Erase all tracking data on which matches have already been output.
	void
	reset() override;

private:
	core::Size n_geometric_constraints_;
	core::Size limit_per_rotamer_combo_;
	RotamerComboCountMap count_per_rotamer_combo_;

};

}
}
}

#endif
