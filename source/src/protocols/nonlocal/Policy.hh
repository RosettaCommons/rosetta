// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/nonlocal/Policy.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_NONLOCAL_POLICY_HH
#define INCLUDED_PROTOCOLS_NONLOCAL_POLICY_HH

// Unit headers
#include <protocols/nonlocal/Policy.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace nonlocal {

/// @class An abstract base class that defines the interface for choosing a
/// fragment from set of possibilities. This class is based on the observation
/// that fragment insertion methods differ primarily in the manner in which they
/// choose from among a set of possibilities.
class Policy : public utility::pointer::ReferenceCount {
	typedef core::fragment::FragSetCOP FragSetCOP;

public:
	/// @brief Provides derived classes with the opportunity to precompute various
	/// properties of the fragment set from which they will have to make choices.
	explicit Policy(FragSetCOP fragments) : fragments_(fragments) {}

	/// @brief Selects a single fragment from a set of possibilities given the
	/// current status of the pose.
	virtual core::Size choose(const core::fragment::Frame& frame,
		const core::pose::Pose& pose) = 0;

	// -- Accessors -- //
	FragSetCOP fragments() const {
		return fragments_;
	}

private:
	FragSetCOP fragments_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_POLICY_HH_
