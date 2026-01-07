// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/DockingPartners.hh
/// @brief  A simple utility class to represent sets of docking chains
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_pose_DockingPartners_hh
#define INCLUDED_core_pose_DockingPartners_hh

#include <utility/vector1.hh>
#include <string>

namespace core {
namespace pose {

//////////////////////////////////////////////////////////////////////////////////
// Docking chain specifications

/// @brief A docking partner specification, based on chain letters
///
struct DockingPartners {
	// No internal invariants, so we're a struct
	utility::vector1<std::string> partner1;
	utility::vector1<std::string> partner2;

	////////////////////////
	/// Utility functions

	/// @brief Split a standard docking partner specification string (e.g "AB_HL") into a DockingPartners specification
	/// @details Right now this assumes single chain letters, but we centralize things here to allow
	/// future format expansion for multi-chain letter.
	static
	DockingPartners
	docking_partners_from_string(std::string const & partner_string);

	/// @brief Return whether the partner specification hasn't been set.
	bool
	is_empty() const;

	/// @brief Return whether the partner specification has both sides filled out.
	bool
	has_both() const;

	/// @brief Return a string of a similar format to that used by docking_partners_from_string()
	/// @details This currently fails silently if there are any multi-chain strings -- don't rely on the round-trip.
	std::string
	str() const;

	/// @brief Return the part of the representational string for partner1
	std::string
	partner1_str() const;

	/// @brief Return the part of the representational string for partner2
	std::string
	partner2_str() const;

	friend
	std::ostream & operator<<( std::ostream & output, DockingPartners const & object_to_output );

	// Comparison for being a map key
	bool
	operator<( DockingPartners const & other ) const;

	bool
	operator==( DockingPartners const & other ) const;
};

} // pose
} // core

#endif // INCLUDED_core_pose_DockingPartners_hh
