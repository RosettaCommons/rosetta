// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/PolicyFactory.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/PolicyFactory.hh>

// C/C++ headers
#include <string>

// External headers
#include <boost/algorithm/string/case_conv.hpp>

// Utility headers
#include <utility/exit.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/util.hh>

// Package headers
#include <protocols/nonlocal/Policy.hh>
#include <protocols/nonlocal/SmoothPolicy.hh>
#include <protocols/nonlocal/UniformPolicy.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace nonlocal {

PolicyOP PolicyFactory::get_policy(const std::string& policy_name,
	core::fragment::FragSetCOP fragments,
	core::Size num_fragments) {
	assert(fragments);
	assert(num_fragments > 0);

	std::string type(policy_name);
	boost::to_lower(type);

	// Operate on a copy of the input to prevent unexpected (and unwanted)
	// modification to the user's fragment data.
	core::fragment::FragSetOP reduced_fragments = fragments->clone();

	// Only consider the top <num_fragments> fragments within each Frame.
	core::fragment::retain_top(num_fragments, reduced_fragments);

	if ( type == "uniform" ) {
		return PolicyOP( new UniformPolicy(reduced_fragments) );
	} else if ( type == "smooth" ) {
		return PolicyOP( new SmoothPolicy(reduced_fragments) );
	} else {
		utility_exit_with_message("Invalid policy_name: " + policy_name);
	}
	return nullptr;
}

}  // namespace nonlocal
}  // namespace protocols
