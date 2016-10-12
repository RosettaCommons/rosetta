// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/ContactCounts.hh
/// @brief  Contact counts.
/// @author David E. Kim (dekim@u.washington.edu)

#ifndef INCLUDED_protocols_frag_picker_ContactCounts_hh
#define INCLUDED_protocols_frag_picker_ContactCounts_hh

// unit headers
#include <protocols/frag_picker/ContactCounts.fwd.hh>

#include <protocols/frag_picker/ContactTypes.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>

#include <string>
#include <map>

// C++ headers

namespace protocols {
namespace frag_picker {

class ContactCounts : public utility::pointer::ReferenceCount {
public:

	ContactCounts() {}

	~ContactCounts() override = default;

	void iterate_neighbor( std::pair<core::Size,core::Size> & query_pair, std::pair<core::Size,core::Size> & neighbor_pair ) {
		neighbor_counts_[query_pair][neighbor_pair]++;
	}

	void iterate( std::pair<core::Size,core::Size> & query_pair ) {
		counts_[query_pair]++;
	}

	std::map<std::pair<core::Size,core::Size>, core::Size> & counts() {
		return counts_;
	}

	bool neighbor_counts_exist( std::pair<core::Size,core::Size> & query_pair ) {
		return (neighbor_counts_.find( query_pair ) != neighbor_counts_.end()) ? true : false;
	}

	std::map<std::pair<core::Size,core::Size>, core::Size> & neighbor_counts( std::pair<core::Size,core::Size> & query_pair ) {
		return neighbor_counts_[query_pair];
	}

private:
	std::map<std::pair<core::Size,core::Size>, core::Size> counts_;
	std::map<std::pair<core::Size,core::Size>, std::map<std::pair<core::Size,core::Size>, core::Size> > neighbor_counts_;
};


} // namespace frag_picker
} // namespace protocols

#endif // INCLUDED_protocols_frag_picker_contact_counts_HH
