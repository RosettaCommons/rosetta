// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/UniformPolicy.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_NONLOCAL_UNIFORMPOLICY_HH
#define INCLUDED_PROTOCOLS_NONLOCAL_UNIFORMPOLICY_HH

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Package headers
#include <protocols/nonlocal/Policy.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace nonlocal {

/// @class Implements the Policy interface. Chooses uniformly among the set of
/// possible fragments at a given position.
class UniformPolicy : public Policy {
	typedef core::fragment::FragSetCOP FragSetCOP;

public:
	/// @brief Provides derived classes with the opportunity to precompute various
	/// properties of the fragment set from which they will have to make choices.
	explicit UniformPolicy(FragSetCOP fragments);

	/// @brief Selects uniformly among the set of possible fragments in <frame>
	virtual core::Size choose(const core::fragment::Frame& frame,
		const core::pose::Pose&);
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_UNIFORMPOLICY_HH_
