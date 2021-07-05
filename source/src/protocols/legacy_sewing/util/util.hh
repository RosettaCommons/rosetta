// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file util.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocol_legacy_sewing_util_util_HH
#define INCLUDED_protocol_legacy_sewing_util_util_HH

//Package headers

//Protocol headers

//Utility

//C++
#include <map>

#include <core/id/AtomID.fwd.hh> // AUTO IWYU For AtomID

namespace protocols {
namespace legacy_sewing  {

std::map<core::id::AtomID, core::id::AtomID>
largest_continuous_atom_map(
	std::map<core::id::AtomID, core::id::AtomID> const & atom_map
);

} //sewing namespace
} //protocols namespace

#endif
