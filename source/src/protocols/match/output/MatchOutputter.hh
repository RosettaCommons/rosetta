// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
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

#ifndef INCLUDED_protocols_match_output_MatchOutputter_hh
#define INCLUDED_protocols_match_output_MatchOutputter_hh

// Unit headers
#include <protocols/match/output/MatchOutputter.fwd.hh>

// Package headers
#include <protocols/match/output/MatchProcessor.hh>

namespace protocols {
namespace match {
namespace output {

class MatchOutputter : public MatchProcessor {
public:
	MatchOutputter();

	virtual
	~MatchOutputter();

	void
	begin_processing();

	void
	end_processing();

	virtual
	void
	process_match(
		match const & m
	);

	virtual
	void
	process_match(
		match_dspos1 const & m
	);

private:

};

}
}
}

#endif
